#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <thread>
#include <filesystem> // For checking file existence
#include <atomic>     // For thread-safe counter
#include <chrono>     // For pauses in progress thread
#include <cstdio>     // for snprintf

// --- Platform-dependent memory utilities ---
#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>
// Note for Windows compilation: need to link with Psapi.lib
// For example, for g++: g++ ... -lpsapi
#elif defined(__linux__)
#include <unistd.h>
#include <ios>
#include <istream>
#elif defined(__APPLE__)
#include <mach/mach.h>
#endif

// --- Utilities ---

// Returns current memory usage (RSS) in megabytes
double getMemoryUsage() {
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS_EX pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc))) {
        return static_cast<double>(pmc.WorkingSetSize) / (1024 * 1024);
    }
    return 0.0;
#elif defined(__linux__)
    std::ifstream statm("/proc/self/statm");
    if (!statm) return 0.0;
    long long rss_pages;
    statm.ignore(std::numeric_limits<std::streamsize>::max(), ' '); // Skip VmSize
    statm >> rss_pages; // Read RSS
    statm.close();
    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024;
    return (static_cast<double>(rss_pages) * page_size_kb) / 1024.0;
#elif defined(__APPLE__)
    mach_task_basic_info_data_t info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &infoCount) != KERN_SUCCESS) {
        return 0.0;
    }
    return static_cast<double>(info.resident_size) / (1024 * 1024);
#else
    return 0.0; // Unsupported platform
#endif
}

// Helper function for formatting time in MM:SS
std::string formatTime(long long seconds) {
    if (seconds < 0) return "--:--";
    long long minutes = seconds / 60;
    seconds %= 60;
    char buffer[32]; // Increased buffer size to prevent overflow
    snprintf(buffer, sizeof(buffer), "%02lld:%02lld", minutes, seconds);
    return buffer;
}


// Function to display progress bar with additional information
void printProgressBar(double percentage, const std::string& message, std::chrono::steady_clock::time_point start_time) {
    percentage = std::max(0.0, std::min(1.0, percentage));
    int val = static_cast<int>(percentage * 100);
    int lpad = static_cast<int>(percentage * 35); // Reduced for additional info
    int rpad = 35 - lpad;

    auto now = std::chrono::steady_clock::now();
    double elapsed_sec = std::chrono::duration_cast<std::chrono::duration<double>>(now - start_time).count();
    
    long long eta_sec = -1;
    if (percentage > 0.001) { // Avoid division by zero and unrealistic ETA at start
        eta_sec = static_cast<long long>(elapsed_sec * (1.0 / percentage - 1.0));
    }

    double mem_usage = getMemoryUsage();
    
    char buffer[256];
    snprintf(buffer, sizeof(buffer), "\r%-20s [%-*s%*s] %3d%% | Elapsed: %s | ETA: %s | Mem: %4.0f MB",
        message.c_str(),
        lpad, std::string(lpad, '=').c_str(),
        rpad, "",
        val,
        formatTime(static_cast<long long>(elapsed_sec)).c_str(),
        formatTime(eta_sec).c_str(),
        mem_usage);

    std::cerr << buffer << std::flush;

    if (percentage >= 1.0) {
        std::cerr << std::endl;
    }
}

// --- Functions for working with K-mers ---

inline int charToBits(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;
    }
}

std::string uint64ToKmer(uint64_t kmer_val, int k) {
    std::string kmer_str(k, ' ');
    uint64_t mask = 3;
    for (int i = 0; i < k; ++i) {
        int bits = (kmer_val >> (2 * (k - 1 - i))) & mask;
        switch (bits) {
            case 0: kmer_str[i] = 'A'; break;
            case 1: kmer_str[i] = 'C'; break;
            case 2: kmer_str[i] = 'G'; break;
            case 3: kmer_str[i] = 'T'; break;
        }
    }
    return kmer_str;
}

// --- Main analyzer class ---

struct ContigInfo {
    std::string name;
    uint64_t length;
    uint64_t global_offset;
};

struct KmerData {
    std::vector<uint64_t> positions; 
    std::vector<uint64_t> distances;
    uint32_t frequency = 0;
};

class KmerAnalyzer {
public:
    KmerAnalyzer(int k_size) : k(k_size) {
        if (k <= 0 || k > 31) {
            throw std::invalid_argument("K-mer size must be between 1 and 31.");
        }
    }

    void process(const std::string& input_fasta, const std::string& output_text_top100, const std::string& output_binary_full, const std::string& contig_index_path, const std::string& cache_path, unsigned int num_threads) {
        bool cache_loaded = loadCacheFromFile(cache_path);

        if (!cache_loaded) {
            readMultiFasta(input_fasta);
            if(genome.empty()){
                std::cerr << "Warning: genome is empty or failed to read file. Analysis stopped." << std::endl;
            } else {
                buildIndexParallel(num_threads);
                calculateAllDistances();
                saveCacheToFile(cache_path);
            }
        }
        
        sortKmersByFrequency();
        writeTopKmers(output_text_top100, 100);
        writeFullBinaryOutput(output_binary_full);
        writeContigIndex(contig_index_path);
    }
    
private:
    const int k;
    std::string genome;
    std::vector<ContigInfo> contig_map;
    std::unordered_map<uint64_t, KmerData> kmer_index;
    std::vector<std::pair<uint64_t, uint32_t>> sorted_kmer_freqs;

    // --- Functions for saving/loading index ---
    const uint32_t MAGIC_NUMBER = 0x4B4D5258; // "KMRX"
    const uint16_t VERSION = 3; // Version 3 with uint64_t and contig information

    void saveCacheToFile(const std::string& filename) {
        auto start_time = std::chrono::steady_clock::now();
        std::ofstream out(filename, std::ios::binary);
        if (!out) { std::cerr << "\nWarning: failed to save cache file: " << filename << std::endl; return; }
        printProgressBar(0.0, "Saving cache", start_time);
        
        out.write(reinterpret_cast<const char*>(&MAGIC_NUMBER), sizeof(MAGIC_NUMBER));
        out.write(reinterpret_cast<const char*>(&VERSION), sizeof(VERSION));
        
        uint64_t contig_count = contig_map.size();
        out.write(reinterpret_cast<const char*>(&contig_count), sizeof(contig_count));
        for(const auto& ci : contig_map) {
            uint32_t name_len = ci.name.length();
            out.write(reinterpret_cast<const char*>(&name_len), sizeof(name_len));
            out.write(ci.name.c_str(), name_len);
            out.write(reinterpret_cast<const char*>(&ci.length), sizeof(ci.length));
            out.write(reinterpret_cast<const char*>(&ci.global_offset), sizeof(ci.global_offset));
        }

        uint64_t map_size = kmer_index.size();
        out.write(reinterpret_cast<const char*>(&map_size), sizeof(map_size));
        uint64_t written_count = 0;
        for (const auto& pair : kmer_index) {
            out.write(reinterpret_cast<const char*>(&pair.first), sizeof(pair.first));
            const auto& data = pair.second;
            out.write(reinterpret_cast<const char*>(&data.frequency), sizeof(data.frequency));
            uint64_t pos_size = data.positions.size();
            out.write(reinterpret_cast<const char*>(&pos_size), sizeof(pos_size));
            out.write(reinterpret_cast<const char*>(data.positions.data()), pos_size * sizeof(uint64_t));
            uint64_t dist_size = data.distances.size();
            out.write(reinterpret_cast<const char*>(&dist_size), sizeof(dist_size));
            out.write(reinterpret_cast<const char*>(data.distances.data()), dist_size * sizeof(uint64_t));
            written_count++;
            if (written_count % 5000 == 0) printProgressBar(static_cast<double>(written_count) / map_size, "Saving cache", start_time);
        }
        printProgressBar(1.0, "Saving cache", start_time);
    }

    bool loadCacheFromFile(const std::string& filename) {
        if (!std::filesystem::exists(filename)) return false;
        std::ifstream in(filename, std::ios::binary);
        if (!in) return false;
        
        auto start_time = std::chrono::steady_clock::now();
        printProgressBar(0.0, "Loading cache", start_time);

        uint32_t magic; uint16_t version;
        in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
        in.read(reinterpret_cast<char*>(&version), sizeof(version));
        if (!in || magic != MAGIC_NUMBER || version != VERSION) { std::cerr << "\nWarning: cache version incompatible." << std::endl; return false; }

        uint64_t contig_count;
        in.read(reinterpret_cast<char*>(&contig_count), sizeof(contig_count));
        contig_map.resize(contig_count);
        for(uint64_t i = 0; i < contig_count; ++i) {
            uint32_t name_len;
            in.read(reinterpret_cast<char*>(&name_len), sizeof(name_len));
            contig_map[i].name.resize(name_len);
            in.read(&contig_map[i].name[0], name_len);
            in.read(reinterpret_cast<char*>(&contig_map[i].length), sizeof(contig_map[i].length));
            in.read(reinterpret_cast<char*>(&contig_map[i].global_offset), sizeof(contig_map[i].global_offset));
        }

        uint64_t map_size;
        in.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
        if (!in) return false;
        if (map_size == 0) { printProgressBar(1.0, "Loading cache", start_time); return true; }
        
        kmer_index.reserve(map_size);
        for (uint64_t i = 0; i < map_size; ++i) {
            uint64_t kmer_val; KmerData data;
            in.read(reinterpret_cast<char*>(&kmer_val), sizeof(kmer_val));
            in.read(reinterpret_cast<char*>(&data.frequency), sizeof(data.frequency));
            uint64_t pos_size, dist_size;
            in.read(reinterpret_cast<char*>(&pos_size), sizeof(pos_size));
            data.positions.resize(pos_size);
            in.read(reinterpret_cast<char*>(data.positions.data()), pos_size * sizeof(uint64_t));
            in.read(reinterpret_cast<char*>(&dist_size), sizeof(dist_size));
            data.distances.resize(dist_size);
            in.read(reinterpret_cast<char*>(data.distances.data()), dist_size * sizeof(uint64_t));
            if (!in) { kmer_index.clear(); return false; }
            kmer_index[kmer_val] = std::move(data);
            if (i % 5000 == 0) printProgressBar(static_cast<double>(i) / map_size, "Loading cache", start_time);
        }
        printProgressBar(1.0, "Loading cache", start_time);
        return true;
    }

    void readMultiFasta(const std::string& filename) {
        auto start_time = std::chrono::steady_clock::now();
        std::ifstream file(filename);
        if (!file) throw std::runtime_error("Cannot open genome file: " + filename);
        
        contig_map.clear();
        genome.clear();
        std::string current_seq_name;
        std::string current_seq;
        uint64_t current_offset = 0;

        auto save_contig = [&]() {
            if (!current_seq_name.empty() && !current_seq.empty()) {
                contig_map.push_back({current_seq_name, static_cast<uint64_t>(current_seq.length()), current_offset});
                genome += current_seq;
                current_offset += current_seq.length();
                current_seq.clear();
            }
        };

        std::string line;
        printProgressBar(0.0, "Reading Multifasta", start_time);
        while (std::getline(file, line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                save_contig();
                size_t first_space = line.find(' ');
                current_seq_name = line.substr(1, first_space - 1);
            } else {
                if (line.back() == '\r') line.pop_back();
                current_seq += line;
            }
        }
        save_contig(); // Save last contig
        printProgressBar(1.0, "Reading Multifasta", start_time);
    }

    void buildIndexParallel(unsigned int num_threads) {
        if (genome.length() < static_cast<size_t>(k)) return;
        auto stage_start_time = std::chrono::steady_clock::now();
        std::string msg = "Indexing (" + std::to_string(num_threads) + " thr.)";
        std::vector<std::thread> threads;
        std::vector<std::unordered_map<uint64_t, std::vector<uint64_t>>> local_indices(num_threads);
        std::atomic<size_t> processed_bases{0};
        std::atomic<bool> building_done{false};
        std::thread progress_thread([&]() {
            while (!building_done) {
                printProgressBar(static_cast<double>(processed_bases) / genome.length(), msg, stage_start_time);
                std::this_thread::sleep_for(std::chrono::milliseconds(200));
            }
        });
        size_t chunk_size = genome.length() / num_threads;
        for (unsigned int i = 0; i < num_threads; ++i) {
            size_t start = i * chunk_size;
            size_t end = (i == num_threads - 1) ? genome.length() : (i + 1) * chunk_size;
            threads.emplace_back(&KmerAnalyzer::buildIndexChunk, this, start, end, std::ref(local_indices[i]), std::ref(processed_bases));
        }
        for (auto& t : threads) t.join();
        building_done = true;
        progress_thread.join();
        printProgressBar(1.0, msg, stage_start_time);
        
        auto reduce_start_time = std::chrono::steady_clock::now();
        printProgressBar(0.0, "Merging maps", reduce_start_time);
        for(size_t i = 0; i < local_indices.size(); ++i) {
            for (const auto& pair : local_indices[i]) {
                kmer_index[pair.first].positions.insert(kmer_index[pair.first].positions.end(), pair.second.begin(), pair.second.end());
            }
            printProgressBar(static_cast<double>(i + 1) / local_indices.size(), "Merging maps", reduce_start_time);
        }
        
        auto final_sort_start_time = std::chrono::steady_clock::now();
        printProgressBar(0.0, "Sorting positions", final_sort_start_time);
        if (!kmer_index.empty()) {
            size_t processed_kmers = 0;
            for(auto& pair : kmer_index) {
                std::sort(pair.second.positions.begin(), pair.second.positions.end());
                pair.second.frequency = pair.second.positions.size();
                processed_kmers++;
                if (processed_kmers % 2000 == 0) printProgressBar(static_cast<double>(processed_kmers) / kmer_index.size(), "Sorting positions", final_sort_start_time);
            }
        }
        printProgressBar(1.0, "Sorting positions", final_sort_start_time);
    }

    void buildIndexChunk(size_t start, size_t end, std::unordered_map<uint64_t, std::vector<uint64_t>>& local_index, std::atomic<size_t>& progress_counter) {
        uint64_t current_kmer = 0;
        uint64_t mask = (1ULL << (2 * k)) - 1;
        int valid_chars_count = 0;
        size_t local_processed_count = 0;
        const size_t update_chunk_size = 4096;
        for (size_t i = start; i < end; ++i) {
            int bits = charToBits(genome[i]);
            if (bits != -1) {
                current_kmer = ((current_kmer << 2) | bits) & mask;
                valid_chars_count++;
            } else {
                valid_chars_count = 0;
            }
            if (valid_chars_count >= k) local_index[current_kmer].push_back(i - k + 1);
            local_processed_count++;
            if(local_processed_count >= update_chunk_size){
                progress_counter.fetch_add(local_processed_count, std::memory_order_relaxed);
                local_processed_count = 0;
            }
        }
        progress_counter.fetch_add(local_processed_count, std::memory_order_relaxed);
    }

    void calculateAllDistances() {
        auto start_time = std::chrono::steady_clock::now();
        if (kmer_index.empty()) { printProgressBar(1.0, "Calculating distances", start_time); return; }
        printProgressBar(0.0, "Calculating distances", start_time);
        size_t processed_count = 0;
        for (auto& pair : kmer_index) {
            auto& data = pair.second;
            if (data.positions.size() > 1) {
                data.distances.reserve(data.positions.size() - 1);
                for (size_t i = 0; i < data.positions.size() - 1; ++i) {
                    data.distances.push_back(data.positions[i + 1] - data.positions[i]);
                }
            }
            processed_count++;
            if (processed_count % 1000 == 0) printProgressBar(static_cast<double>(processed_count) / kmer_index.size(), "Calculating distances", start_time);
        }
        printProgressBar(1.0, "Calculating distances", start_time);
    }

    void sortKmersByFrequency() {
        auto start_time = std::chrono::steady_clock::now();
        printProgressBar(0.0, "Sorting k-mers", start_time);
        for (const auto& pair : kmer_index) sorted_kmer_freqs.push_back({pair.first, pair.second.frequency});
        std::sort(sorted_kmer_freqs.begin(), sorted_kmer_freqs.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
        printProgressBar(1.0, "Sorting k-mers", start_time);
    }
    
    void writeTopKmers(const std::string& filename, int top_n) {
        auto start_time = std::chrono::steady_clock::now();
        std::ofstream out(filename);
        if (!out) throw std::runtime_error("Cannot open file for writing: " + filename);
        out << "# Top " << top_n << " most frequent " << k << "-mers\n";
        out << "# Format: KMER\tFREQUENCY\tDISTANCES(pos2-pos1,pos3-pos2,...)\n";
        int count_to_write = std::min(top_n, static_cast<int>(sorted_kmer_freqs.size()));
        if (count_to_write == 0) { printProgressBar(1.0, "Writing top-100 (txt)", start_time); return; }
        printProgressBar(0.0, "Writing top-100 (txt)", start_time);
        for (int i = 0; i < count_to_write; ++i) {
            uint64_t kmer_val = sorted_kmer_freqs[i].first;
            const auto& data = kmer_index.at(kmer_val);
            out << uint64ToKmer(kmer_val, k) << "\t" << data.frequency << "\t";
            if (!data.distances.empty()) {
                for (size_t j = 0; j < data.distances.size(); ++j) out << data.distances[j] << (j == data.distances.size() - 1 ? "" : ",");
            }
            out << "\n";
            printProgressBar(static_cast<double>(i+1)/count_to_write, "Writing top-100 (txt)", start_time);
        }
        printProgressBar(1.0, "Writing top-100 (txt)", start_time);
    }

    void writeFullBinaryOutput(const std::string& filename) {
        auto start_time = std::chrono::steady_clock::now();
        std::ofstream out(filename, std::ios::binary);
        if (!out) throw std::runtime_error("Cannot open full binary file for writing: " + filename);
        
        printProgressBar(0.0, "Writing (bin-full)", start_time);
        uint64_t map_size = sorted_kmer_freqs.size();
        out.write(reinterpret_cast<const char*>(&map_size), sizeof(map_size));
        uint64_t written_count = 0;
        for (const auto& sorted_pair : sorted_kmer_freqs) {
            uint64_t kmer_val = sorted_pair.first;
            const auto& data = kmer_index.at(kmer_val);
            out.write(reinterpret_cast<const char*>(&kmer_val), sizeof(kmer_val));
            out.write(reinterpret_cast<const char*>(&data.frequency), sizeof(data.frequency));
            uint64_t pos_size = data.positions.size();
            out.write(reinterpret_cast<const char*>(&pos_size), sizeof(pos_size));
            out.write(reinterpret_cast<const char*>(data.positions.data()), pos_size * sizeof(uint64_t));
            uint64_t dist_size = data.distances.size();
            out.write(reinterpret_cast<const char*>(&dist_size), sizeof(dist_size));
            out.write(reinterpret_cast<const char*>(data.distances.data()), dist_size * sizeof(uint64_t));
            written_count++;
            if (written_count % 2000 == 0) printProgressBar(static_cast<double>(written_count) / map_size, "Writing (bin-full)", start_time);
        }
        printProgressBar(1.0, "Writing (bin-full)", start_time);
    }

    void writeContigIndex(const std::string& filename) {
        auto start_time = std::chrono::steady_clock::now();
        std::ofstream out(filename);
        if (!out) { std::cerr << "\nWarning: failed to write contig index: " << filename << std::endl; return; }
        printProgressBar(0.0, "Writing .cidx file", start_time);
        out << "#contig_name\tlength\tglobal_offset\n";
        for(const auto& ci : contig_map) {
            out << ci.name << "\t" << ci.length << "\t" << ci.global_offset << "\n";
        }
        printProgressBar(1.0, "Writing .cidx file", start_time);
    }
};

int main(int argc, char* argv[]) {
    if (argc < 5 || argc > 7) {
        std::cerr << "Usage: " << argv[0] << " <in.fasta> <out_top100.txt> <out_full.bin> <out_contigs.cidx> [cache.kidx] [threads]" << std::endl;
        std::cerr << "  <out_contigs.cidx> - text file with contig information." << std::endl;
        std::cerr << "  [cache.kidx]       - optional. Cache file for saving/loading." << std::endl;
        std::cerr << "  [num_threads]      - optional. Number of threads (default: all available)." << std::endl;
        return 1;
    }

    std::string input_file = argv[1];
    std::string output_text_file = argv[2];
    std::string output_binary_file = argv[3];
    std::string contig_index_file = argv[4];
    std::string cache_file = (argc >= 6) ? argv[5] : "kmer_cache.kidx";
    
    unsigned int num_threads = std::thread::hardware_concurrency();
    if (argc >= 7) {
        try {
            num_threads = std::stoul(argv[6]);
            if (num_threads == 0) num_threads = 1;
        } catch (...) { std::cerr << "Error: invalid value for num_threads." << std::endl; return 1; }
    }
    
    const int K_SIZE = 13;

    try {
        KmerAnalyzer analyzer(K_SIZE);
        analyzer.process(input_file, output_text_file, output_binary_file, contig_index_file, cache_file, num_threads);
        std::cout << "\nAnalysis completed successfully." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

