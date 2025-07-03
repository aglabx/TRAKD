#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <chrono>
#include <iomanip> // for std::setprecision
#include <sstream> // for std::stringstream
#include <thread>
#include <atomic>
#include <mutex>

// --- Platform-dependent memory utilities ---
#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>
#elif defined(__linux__)
#include <unistd.h>
#include <ios>
#include <istream>
#elif defined(__APPLE__)
#include <mach/mach.h>
#endif

// --- Utilities ---

double getMemoryUsage() {
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS_EX pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc))) return static_cast<double>(pmc.WorkingSetSize) / (1024 * 1024);
    return 0.0;
#elif defined(__linux__)
    std::ifstream statm("/proc/self/statm");
    if (!statm) return 0.0;
    long long rss_pages;
    statm.ignore(std::numeric_limits<std::streamsize>::max(), ' ');
    statm >> rss_pages;
    statm.close();
    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024;
    return (static_cast<double>(rss_pages) * page_size_kb) / 1024.0;
#elif defined(__APPLE__)
    mach_task_basic_info_data_t info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &infoCount) != KERN_SUCCESS) return 0.0;
    return static_cast<double>(info.resident_size) / (1024 * 1024);
#else
    return 0.0;
#endif
}

std::string formatTime(long long seconds) {
    if (seconds < 0) return "--:--";
    long long minutes = seconds / 60;
    seconds %= 60;
    char buffer[32];
    snprintf(buffer, sizeof(buffer), "%02lld:%02lld", minutes, seconds);
    return buffer;
}

void printProgressBar(double percentage, const std::string& message, std::chrono::steady_clock::time_point start_time, const std::string& extra_info = "") {
    percentage = std::max(0.0, std::min(1.0, percentage));
    int val = static_cast<int>(percentage * 100);
    int lpad = static_cast<int>(percentage * 40);
    int rpad = 40 - lpad;
    auto now = std::chrono::steady_clock::now();
    double elapsed_sec = std::chrono::duration_cast<std::chrono::duration<double>>(now - start_time).count();
    long long eta_sec = -1;
    if (percentage > 0.001) eta_sec = static_cast<long long>(elapsed_sec * (1.0 / percentage - 1.0));
    double mem_usage = getMemoryUsage();
    char buffer[256];
    snprintf(buffer, sizeof(buffer), "\r%-20s [%-*s%*s] %3d%% %-15s | Elapsed: %s | ETA: %s | Mem: %.0f MB",
        message.c_str(), lpad, std::string(lpad, '=').c_str(), rpad, "", val,
        extra_info.c_str(), formatTime(static_cast<long long>(elapsed_sec)).c_str(), formatTime(eta_sec).c_str(), mem_usage);
    std::cerr << buffer << std::flush;
    if (percentage >= 1.0) std::cerr << std::endl;
}

// --- Main structures ---

struct ContigInfo {
    std::string name;
    uint64_t length;
    uint64_t global_offset;
};

struct KmerRawData {
    uint32_t tf;
    std::vector<uint64_t> positions;
};

struct Locus {
    uint64_t start;
    uint64_t end;
    bool operator<(const Locus& other) const { return start < other.start; }
};

class LocusBedGenerator {
public:
    void generate(const std::string& input_binary_full, const std::string& input_fai, const std::string& output_bed, 
                  uint32_t min_kmer_tf, uint32_t locus_gap_threshold, uint32_t min_locus_kmer_count, unsigned int num_threads) {
        
        loadFaiIndex(input_fai);
        
        std::vector<KmerRawData> kmer_jobs;
        readAllKmerData(input_binary_full, kmer_jobs, min_kmer_tf);
        
        if (kmer_jobs.empty()) {
            std::cerr << "No k-mers for analysis after applying filter." << std::endl;
            std::ofstream out(output_bed);
            out << "track name=\"MergedKmerLoci\" description=\"Merged loci from k-mer chains\"\n";
            return;
        }

        std::vector<Locus> all_loci = findPrimaryLociParallel(kmer_jobs, locus_gap_threshold, min_locus_kmer_count, num_threads);
        std::vector<Locus> merged_loci = mergeLoci(all_loci);
        writeBedFile(output_bed, merged_loci);
    }

private:
    std::vector<ContigInfo> contig_map;

    void loadFaiIndex(const std::string& filename) {
        auto start_time = std::chrono::steady_clock::now();
        printProgressBar(0.0, "Reading .fai index", start_time);
        std::ifstream file(filename);
        if (!file) throw std::runtime_error("Cannot open .fai index file: " + filename);
        
        contig_map.clear();
        std::string line;
        uint64_t current_offset = 0;
        while(std::getline(file, line)) {
            if (line.empty()) continue;
            std::stringstream ss(line);
            ContigInfo info;
            ss >> info.name >> info.length;
            info.global_offset = current_offset;
            contig_map.push_back(info);
            current_offset += info.length;
        }
        if (contig_map.empty()) throw std::runtime_error(".fai index file is empty or has invalid format.");
        printProgressBar(1.0, "Reading .fai index", start_time);
    }

    void readAllKmerData(const std::string& filename, std::vector<KmerRawData>& kmer_jobs, uint32_t min_kmer_tf) {
        auto start_time = std::chrono::steady_clock::now();
        std::ifstream in(filename, std::ios::binary);
        if (!in) throw std::runtime_error("Cannot open full binary file: " + filename);
        
        uint64_t total_kmers_in_file;
        in.read(reinterpret_cast<char*>(&total_kmers_in_file), sizeof(total_kmers_in_file));
        if (!in || total_kmers_in_file == 0) return;
        
        std::cerr << "Header information: " << total_kmers_in_file << " k-mers in file." << std::endl;
        printProgressBar(0.0, "Reading data", start_time);

        for (uint64_t i = 0; i < total_kmers_in_file; ++i) {
            uint64_t kmer_val;
            uint32_t kmer_tf;
            in.read(reinterpret_cast<char*>(&kmer_val), sizeof(kmer_val));
            in.read(reinterpret_cast<char*>(&kmer_tf), sizeof(kmer_tf));

            if (kmer_tf < min_kmer_tf) {
                std::cerr << "\nReached k-mer frequency threshold (" << min_kmer_tf << "). Loaded " << i << " k-mers." << std::endl;
                break;
            }

            uint64_t pos_size, dist_size;
            in.read(reinterpret_cast<char*>(&pos_size), sizeof(pos_size));
            
            // *** DETAILED LOGGING AND SAFETY CHECK ***
            if (i % 1000 == 0) { // Log not for every k-mer to avoid cluttering output
                std::cerr << "\n[DEBUG] Reading k-mer " << i << ": tf=" << kmer_tf << ", pos_size=" << pos_size;
            }
            const uint64_t SANE_LIMIT = 2000000000; // 2 billion positions - very generous limit
            if (pos_size > SANE_LIMIT) {
                std::cerr << "\n\nERROR: Found k-mer with unrealistically large number of positions (" << pos_size << ").\n"
                          << "File " << filename << " may be corrupted or have incorrect format.\n"
                          << "Emergency exit to avoid memory allocation failure." << std::endl;
                exit(1);
            }
            
            KmerRawData job;
            job.tf = kmer_tf;
            try {
                job.positions.resize(pos_size);
            } catch (const std::bad_alloc& e) {
                 std::cerr << "\n\nCRITICAL ERROR: std::bad_alloc when trying to allocate memory for " << pos_size << " positions.\n"
                           << "k-mer index: " << i << ", tf: " << kmer_tf << std::endl;
                 throw;
            }
            
            in.read(reinterpret_cast<char*>(job.positions.data()), pos_size * sizeof(uint64_t));
            in.read(reinterpret_cast<char*>(&dist_size), sizeof(dist_size));
            in.seekg(dist_size * sizeof(uint64_t), std::ios_base::cur);
            
            if (!in) throw std::runtime_error("Error reading binary file.");
            kmer_jobs.push_back(std::move(job));

            if (i % 500 == 0) {
                printProgressBar(static_cast<double>(i + 1) / total_kmers_in_file, "Reading data", start_time);
            }
        }
        printProgressBar(1.0, "Reading data", start_time);
    }

    void findLociWorker(const std::vector<KmerRawData>& jobs, std::vector<Locus>& local_loci, std::atomic<size_t>& next_job_idx, uint32_t locus_gap_threshold, uint32_t min_locus_kmer_count) {
        while (true) {
            size_t i = next_job_idx.fetch_add(1);
            if (i >= jobs.size()) break;
            const auto& job = jobs[i];
            if (job.positions.empty()) continue;
            uint64_t current_locus_start = job.positions[0];
            uint32_t current_locus_kmer_count = 1;
            size_t contig_idx = findContigIndex(current_locus_start);
            for (size_t p_idx = 1; p_idx < job.positions.size(); ++p_idx) {
                uint64_t current_pos = job.positions[p_idx];
                uint64_t prev_pos = job.positions[p_idx - 1];
                bool crosses_boundary = (findContigIndex(current_pos) != contig_idx);
                if (crosses_boundary || (current_pos - prev_pos > locus_gap_threshold)) {
                    if (current_locus_kmer_count >= min_locus_kmer_count) local_loci.push_back({current_locus_start, prev_pos});
                    current_locus_start = current_pos;
                    current_locus_kmer_count = 1;
                    contig_idx = findContigIndex(current_pos);
                } else {
                    current_locus_kmer_count++;
                }
            }
            if (current_locus_kmer_count >= min_locus_kmer_count) local_loci.push_back({current_locus_start, job.positions.back()});
        }
    }

    std::vector<Locus> findPrimaryLociParallel(const std::vector<KmerRawData>& kmer_jobs, uint32_t locus_gap_threshold, uint32_t min_locus_kmer_count, unsigned int num_threads) {
        auto start_time = std::chrono::steady_clock::now();
        std::string msg = "Finding loci (" + std::to_string(num_threads) + " thr.)";
        std::vector<std::thread> threads;
        std::vector<std::vector<Locus>> thread_local_loci(num_threads);
        std::atomic<size_t> next_job_idx{0};
        for (unsigned int i = 0; i < num_threads; ++i) {
            threads.emplace_back(&LocusBedGenerator::findLociWorker, this, std::ref(kmer_jobs), std::ref(thread_local_loci[i]), std::ref(next_job_idx), locus_gap_threshold, min_locus_kmer_count);
        }
        while(true) {
            size_t processed_count = next_job_idx.load();
            if (processed_count > kmer_jobs.size()) processed_count = kmer_jobs.size();
            std::string extra_info = "[" + std::to_string(processed_count) + "/" + std::to_string(kmer_jobs.size()) + "]";
            printProgressBar(static_cast<double>(processed_count) / kmer_jobs.size(), msg, start_time, extra_info);
            if (processed_count >= kmer_jobs.size()) break;
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        for (auto& t : threads) t.join();
        std::vector<Locus> all_loci;
        for (const auto& local_loci : thread_local_loci) all_loci.insert(all_loci.end(), local_loci.begin(), local_loci.end());
        return all_loci;
    }

    std::vector<Locus> mergeLoci(std::vector<Locus>& loci) {
        auto start_time = std::chrono::steady_clock::now();
        if (loci.empty()) { printProgressBar(1.0, "Merging loci", start_time); return {}; }
        printProgressBar(0.0, "Merging loci", start_time, "[1/2 Sorting]");
        std::sort(loci.begin(), loci.end());
        printProgressBar(0.5, "Merging loci", start_time, "[2/2 Merging]");
        std::vector<Locus> merged_loci;
        merged_loci.push_back(loci[0]);
        for (size_t i = 1; i < loci.size(); ++i) {
            if (loci[i].start <= merged_loci.back().end) {
                merged_loci.back().end = std::max(merged_loci.back().end, loci[i].end);
            } else {
                merged_loci.push_back(loci[i]);
            }
            if (i % 10000 == 0) printProgressBar(0.5 + 0.5 * (static_cast<double>(i) / loci.size()), "Merging loci", start_time, "[2/2 Merging]");
        }
        printProgressBar(1.0, "Merging loci", start_time);
        return merged_loci;
    }

    void writeBedFile(const std::string& filename, const std::vector<Locus>& merged_loci) {
        auto start_time = std::chrono::steady_clock::now();
        std::ofstream out(filename);
        if (!out) throw std::runtime_error("Cannot open BED file for writing: " + filename);
        out << "track name=\"MergedKmerLoci\" description=\"Merged loci from k-mer chains\"\n";
        if (merged_loci.empty()) { printProgressBar(1.0, "Writing BED file", start_time); return; }
        printProgressBar(0.0, "Writing BED file", start_time);
        for (size_t i = 0; i < merged_loci.size(); ++i) {
            const auto& locus = merged_loci[i];
            size_t contig_idx = findContigIndex(locus.start);
            const auto& contig = contig_map[contig_idx];
            uint64_t local_start = locus.start - contig.global_offset;
            uint64_t local_end = locus.end - contig.global_offset + 13; // +k for full length
            if (local_end > contig.length) local_end = contig.length;
            out << contig.name << "\t" << local_start << "\t" << local_end << "\t"
                << "locus_" << i+1 << "_len_" << (local_end - local_start) << "\t" << "0\t.\n";
            if (i % 1000 == 0) printProgressBar(static_cast<double>(i+1)/merged_loci.size(), "Writing BED file", start_time);
        }
        printProgressBar(1.0, "Writing BED file", start_time);
    }

    size_t findContigIndex(uint64_t global_pos) const {
        auto it = std::upper_bound(contig_map.begin(), contig_map.end(), global_pos, 
            [](uint64_t pos, const ContigInfo& info){ return pos < info.global_offset; });
        if (it == contig_map.begin()) return 0;
        return std::distance(contig_map.begin(), it) - 1;
    }
};

int main(int argc, char* argv[]) {
    if (argc < 4 || argc > 8) {
        std::cerr << "Usage: " << argv[0] << " <input.bin> <input.fai> <output.bed> [min_tf] [locus_gap] [min_locus_kmers] [threads]" << std::endl;
        std::cerr << "  <input.bin>         - binary file from kmer_analyzer." << std::endl;
        std::cerr << "  <input.fai>         - standard FASTA index (.fai)." << std::endl;
        std::cerr << "  <output.bed>        - BED file with merged loci." << std::endl;
        std::cerr << "  [min_kmer_tf]       - optional. Min. k-mer frequency for analysis (default: 100)." << std::endl;
        std::cerr << "  [locus_gap]         - optional. Max. gap within locus (default: 10000)." << std::endl;
        std::cerr << "  [min_locus_kmers]   - optional. Min. number of k-mers in locus (default: 2)." << std::endl;
        std::cerr << "  [threads]           - optional. Number of threads (default: all available)." << std::endl;
        return 1;
    }

    std::string input_file = argv[1];
    std::string fai_file = argv[2];
    std::string output_file = argv[3];
    uint32_t min_kmer_tf = 100;
    uint32_t locus_gap_threshold = 10000;
    uint32_t min_locus_kmer_count = 2;
    unsigned int num_threads = std::thread::hardware_concurrency();

    if (argc >= 5) { try { min_kmer_tf = std::stoul(argv[4]); } catch (...) { std::cerr << "Error: invalid value for min_kmer_tf." << std::endl; return 1; } }
    if (argc >= 6) { try { locus_gap_threshold = std::stoul(argv[5]); } catch (...) { std::cerr << "Error: invalid value for locus_gap_threshold." << std::endl; return 1; } }
    if (argc >= 7) { try { min_locus_kmer_count = std::stoul(argv[6]); } catch (...) { std::cerr << "Error: invalid value for min_locus_kmers." << std::endl; return 1; } }
    if (argc >= 8) {
        try {
            num_threads = std::stoul(argv[7]);
            if (num_threads == 0) num_threads = 1;
        } catch (...) { std::cerr << "Error: invalid value for num_threads." << std::endl; return 1; }
    }
    
    try {
        LocusBedGenerator merger;
        merger.generate(input_file, fai_file, output_file, min_kmer_tf, locus_gap_threshold, min_locus_kmer_count, num_threads);
        std::cout << "\nBED file creation completed successfully." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

