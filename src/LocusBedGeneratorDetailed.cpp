#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <chrono>
#include <iomanip> // for std::setprecision
#include <sstream> // for std::stringstream
#include <thread>
#include <atomic>
#include <mutex>
#include <cmath> // for log2
#include <numeric> // for std::iota

// --- Utilities ---

// Helper function for formatting time in MM:SS
std::string formatTime(long long seconds) {
    if (seconds < 0) return "--:--";
    long long minutes = seconds / 60;
    seconds %= 60;
    char buffer[32];
    snprintf(buffer, sizeof(buffer), "%02lld:%02lld", minutes, seconds);
    return buffer;
}

// Function to display progress bar
void printProgressBar(double percentage, const std::string& message, std::chrono::steady_clock::time_point start_time, const std::string& extra_info = "") {
    percentage = std::max(0.0, std::min(1.0, percentage));
    int val = static_cast<int>(percentage * 100);
    int lpad = static_cast<int>(percentage * 40);
    int rpad = 40 - lpad;

    auto now = std::chrono::steady_clock::now();
    double elapsed_sec = std::chrono::duration_cast<std::chrono::duration<double>>(now - start_time).count();
    
    long long eta_sec = -1;
    if (percentage > 0.001) {
        eta_sec = static_cast<long long>(elapsed_sec * (1.0 / percentage - 1.0));
    }
    
    char buffer[256];
    snprintf(buffer, sizeof(buffer), "\r%-20s [%-*s%*s] %3d%% %-15s | Elapsed: %s | ETA: %s",
        message.c_str(),
        lpad, std::string(lpad, '=').c_str(),
        rpad, "",
        val,
        extra_info.c_str(),
        formatTime(static_cast<long long>(elapsed_sec)).c_str(),
        formatTime(eta_sec).c_str());

    std::cerr << buffer << std::flush;

    if (percentage >= 1.0) {
        std::cerr << std::endl;
    }
}

// Function to convert 64-bit number to k-mer string
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


// --- Main structures ---

struct ContigInfo {
    std::string name;
    uint64_t length;
    uint64_t global_offset;
};

struct KmerRawData {
    uint64_t kmer_val;
    uint32_t tf;
    std::vector<uint64_t> positions;
    std::vector<uint64_t> distances;
};

struct PrimaryLocus {
    uint64_t start;
    uint64_t end;
    size_t kmer_db_index; 
    bool operator<(const PrimaryLocus& other) const { return start < other.start; }
};

struct MergedLocus {
    uint64_t start;
    uint64_t end;
    std::unordered_set<size_t> contributing_kmer_indices;
    std::vector<double> jaccard_values;
};

struct BedRecord {
    std::string chrom;
    uint64_t start;
    uint64_t end;
    std::string name;
    uint64_t score;
    char strand = '.';
    
    bool operator<(const BedRecord& other) const {
        if (chrom != other.chrom) {
            return chrom < other.chrom;
        }
        return start < other.start;
    }
};


class LocusBedGeneratorDetailed {
public:
    void generate(const std::string& input_binary_full, const std::string& input_fai, const std::string& output_bed, 
                  uint32_t min_kmer_tf, uint32_t locus_gap_threshold, uint32_t min_locus_kmer_count, 
                  double intersection_overlap, double jaccard_threshold, uint32_t min_locus_length, unsigned int num_threads) {
        
        loadFaiIndex(input_fai);
        
        std::vector<KmerRawData> kmer_db;
        readAllKmerData(input_binary_full, kmer_db, min_kmer_tf);
        
        if (kmer_db.empty()) {
            std::cerr << "No k-mers for analysis after applying filter." << std::endl;
            std::ofstream out(output_bed);
            out << "track name=\"DetailedLoci\" description=\"Merged loci with signature k-mer stats\"\n";
            return;
        }

        std::vector<PrimaryLocus> all_primary_loci = findPrimaryLociParallel(kmer_db, locus_gap_threshold, min_locus_kmer_count, num_threads);

        std::vector<MergedLocus> merged_loci = mergeLoci(all_primary_loci, kmer_db, intersection_overlap, jaccard_threshold);
        
        writeBedFile(output_bed, merged_loci, kmer_db, locus_gap_threshold, min_locus_length);
    }

private:
    std::vector<ContigInfo> contig_map;
    const int K_SIZE = 13; 

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

    void readAllKmerData(const std::string& filename, std::vector<KmerRawData>& kmer_db, uint32_t min_kmer_tf) {
        auto start_time = std::chrono::steady_clock::now();
        std::ifstream in(filename, std::ios::binary);
        if (!in) throw std::runtime_error("Cannot open full binary file: " + filename);
        
        uint64_t total_kmers_in_file;
        in.read(reinterpret_cast<char*>(&total_kmers_in_file), sizeof(total_kmers_in_file));
        if (!in || total_kmers_in_file == 0) return;
        
        printProgressBar(0.0, "Reading data", start_time);
        for (uint64_t i = 0; i < total_kmers_in_file; ++i) {
            KmerRawData job;
            in.read(reinterpret_cast<char*>(&job.kmer_val), sizeof(job.kmer_val));
            in.read(reinterpret_cast<char*>(&job.tf), sizeof(job.tf));
            if (job.tf < min_kmer_tf) break;
            uint64_t pos_size, dist_size;
            in.read(reinterpret_cast<char*>(&pos_size), sizeof(pos_size));
            job.positions.resize(pos_size);
            in.read(reinterpret_cast<char*>(job.positions.data()), pos_size * sizeof(uint64_t));
            in.read(reinterpret_cast<char*>(&dist_size), sizeof(dist_size));
            job.distances.resize(dist_size);
            in.read(reinterpret_cast<char*>(job.distances.data()), dist_size * sizeof(uint64_t));
            if (!in) throw std::runtime_error("Error reading binary file.");
            kmer_db.push_back(std::move(job));
            if (i % 500 == 0) printProgressBar(static_cast<double>(i + 1) / total_kmers_in_file, "Reading data", start_time);
        }
        printProgressBar(1.0, "Reading data", start_time);
    }

    void findLociWorker(const std::vector<KmerRawData>& jobs, std::vector<PrimaryLocus>& local_loci, std::atomic<size_t>& next_job_idx, uint32_t locus_gap_threshold, uint32_t min_locus_kmer_count) {
        while (true) {
            size_t i = next_job_idx.fetch_add(1);
            if (i >= jobs.size()) break;
            const auto& job = jobs[i];
            if (job.positions.size() < min_locus_kmer_count) continue;

            uint64_t current_locus_start = job.positions[0];
            uint32_t current_locus_kmer_count = 1;
            size_t contig_idx = findContigIndex(current_locus_start);

            for (size_t p_idx = 1; p_idx < job.positions.size(); ++p_idx) {
                uint64_t current_pos = job.positions[p_idx];
                uint64_t prev_pos = job.positions[p_idx - 1];
                bool crosses_boundary = (findContigIndex(current_pos) != contig_idx);
                if (crosses_boundary || (current_pos - prev_pos > locus_gap_threshold)) {
                    if (current_locus_kmer_count >= min_locus_kmer_count) {
                        local_loci.push_back({current_locus_start, prev_pos, i});
                    }
                    current_locus_start = current_pos;
                    current_locus_kmer_count = 1;
                    contig_idx = findContigIndex(current_pos);
                } else {
                    current_locus_kmer_count++;
                }
            }
            if (current_locus_kmer_count >= min_locus_kmer_count) {
                local_loci.push_back({current_locus_start, job.positions.back(), i});
            }
        }
    }

    std::vector<PrimaryLocus> findPrimaryLociParallel(const std::vector<KmerRawData>& kmer_jobs, uint32_t locus_gap_threshold, uint32_t min_locus_kmer_count, unsigned int num_threads) {
        auto start_time = std::chrono::steady_clock::now();
        std::string msg = "Finding loci (" + std::to_string(num_threads) + " thr.)";
        std::vector<std::thread> threads;
        std::vector<std::vector<PrimaryLocus>> thread_local_loci(num_threads);
        std::atomic<size_t> next_job_idx{0};

        for (unsigned int i = 0; i < num_threads; ++i) {
            threads.emplace_back(&LocusBedGeneratorDetailed::findLociWorker, this, std::ref(kmer_jobs), std::ref(thread_local_loci[i]), std::ref(next_job_idx), locus_gap_threshold, min_locus_kmer_count);
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

        std::vector<PrimaryLocus> all_loci;
        for (const auto& local_loci : thread_local_loci) {
            all_loci.insert(all_loci.end(), local_loci.begin(), local_loci.end());
        }
        return all_loci;
    }
    
    void getPositionalSignature(const std::unordered_set<size_t>& kmer_indices, const std::vector<KmerRawData>& kmer_db, uint64_t start, uint64_t end, std::unordered_set<uint64_t>& signature) {
        signature.clear();
        for (size_t idx : kmer_indices) {
            const auto& positions = kmer_db[idx].positions;
            auto lower = std::lower_bound(positions.begin(), positions.end(), start);
            auto upper = std::upper_bound(positions.begin(), positions.end(), end);
            for (auto it = lower; it != upper; ++it) {
                for(int offset = 0; offset < K_SIZE; ++offset) {
                    signature.insert(*it + offset);
                }
            }
        }
    }

    double calculatePositionalJaccard(const std::unordered_set<size_t>& set1_indices, const std::unordered_set<size_t>& set2_indices, const std::vector<KmerRawData>& kmer_db, uint64_t locus1_start, uint64_t locus1_end, uint64_t locus2_start, uint64_t locus2_end) {
        std::unordered_set<uint64_t> signature1, signature2;
        getPositionalSignature(set1_indices, kmer_db, locus1_start, locus1_end, signature1);
        getPositionalSignature(set2_indices, kmer_db, locus2_start, locus2_end, signature2);

        if (signature1.empty() || signature2.empty()) return 0.0;
        size_t intersection_size = 0;
        const auto& smaller_set = (signature1.size() < signature2.size()) ? signature1 : signature2;
        const auto& larger_set = (signature1.size() < signature2.size()) ? signature2 : signature1;
        for (const auto& pos : smaller_set) {
            if (larger_set.count(pos)) intersection_size++;
        }
        size_t union_size = signature1.size() + signature2.size() - intersection_size;
        return (union_size == 0) ? 0.0 : static_cast<double>(intersection_size) / union_size;
    }

    std::vector<MergedLocus> mergeLoci(std::vector<PrimaryLocus>& primary_loci, const std::vector<KmerRawData>& kmer_db, double intersection_overlap, double jaccard_threshold) {
        auto start_time = std::chrono::steady_clock::now();
        if (primary_loci.empty()) { printProgressBar(1.0, "Merging loci", start_time); return {}; }
        
        printProgressBar(0.0, "Merging loci", start_time, "[1/2 Sorting]");
        std::sort(primary_loci.begin(), primary_loci.end());
        
        printProgressBar(0.5, "Merging loci", start_time, "[2/2 Merging]");
        std::vector<MergedLocus> merged_loci;
        if (primary_loci.empty()) return merged_loci;

        merged_loci.push_back({primary_loci[0].start, primary_loci[0].end, {primary_loci[0].kmer_db_index}, {}});

        for (size_t i = 1; i < primary_loci.size(); ++i) {
            auto& current_merged = merged_loci.back();
            const auto& next_primary = primary_loci[i];

            if (findContigIndex(next_primary.start) != findContigIndex(current_merged.start)) {
                 merged_loci.push_back({next_primary.start, next_primary.end, {next_primary.kmer_db_index}, {}});
                 continue;
            }

            if (next_primary.start <= current_merged.end) {
                uint64_t overlap_start = std::max(current_merged.start, next_primary.start);
                uint64_t overlap_end = std::min(current_merged.end, next_primary.end);
                uint64_t overlap_len = (overlap_end > overlap_start) ? (overlap_end - overlap_start) : 0;
                uint64_t min_locus_len = std::min(current_merged.end - current_merged.start, next_primary.end - next_primary.start);
                double overlap_frac = (min_locus_len > 0) ? static_cast<double>(overlap_len) / min_locus_len : 0.0;

                bool should_merge = false;
                double j = 0.0;
                if (overlap_frac >= intersection_overlap) {
                    should_merge = true;
                    j = 1.0; // Forced merge
                } else {
                    std::unordered_set<size_t> next_set = {next_primary.kmer_db_index};
                    j = calculatePositionalJaccard(current_merged.contributing_kmer_indices, next_set, kmer_db, current_merged.start, current_merged.end, next_primary.start, next_primary.end);
                    if (j >= jaccard_threshold) {
                        should_merge = true;
                    }
                }
                
                if (should_merge) {
                    current_merged.end = std::max(current_merged.end, next_primary.end);
                    current_merged.contributing_kmer_indices.insert(next_primary.kmer_db_index);
                    current_merged.jaccard_values.push_back(j);
                } else {
                    merged_loci.push_back({next_primary.start, next_primary.end, {next_primary.kmer_db_index}, {}});
                }
            } else {
                merged_loci.push_back({next_primary.start, next_primary.end, {next_primary.kmer_db_index}, {}});
            }

            if (i % 10000 == 0) {
                printProgressBar(0.5 + 0.5 * (static_cast<double>(i) / primary_loci.size()), "Merging loci", start_time, "[2/2 Merging]");
            }
        }
        printProgressBar(1.0, "Merging loci", start_time);
        return merged_loci;
    }


    void writeBedFile(const std::string& filename, std::vector<MergedLocus>& merged_loci, const std::vector<KmerRawData>& kmer_db, uint32_t locus_gap_threshold, uint32_t min_locus_length) {
        auto start_time = std::chrono::steady_clock::now();
        std::ofstream out(filename);
        if (!out) throw std::runtime_error("Cannot open BED file for writing: " + filename);
        
        out << "track name=\"DetailedLoci\" description=\"Merged loci with signature k-mer stats\"\n";
        if (merged_loci.empty()) { printProgressBar(1.0, "Writing BED file", start_time); return; }

        std::vector<BedRecord> bed_records;
        
        printProgressBar(0.0, "Preparing to write", start_time);
        for (size_t i = 0; i < merged_loci.size(); ++i) {
            auto& locus = merged_loci[i];
            
            size_t signature_kmer_idx = 0;
            uint32_t max_kmers_in_locus = 0;

            if (!locus.contributing_kmer_indices.empty()){
                signature_kmer_idx = *locus.contributing_kmer_indices.begin();
                for (size_t kmer_idx : locus.contributing_kmer_indices) {
                    const auto& positions = kmer_db[kmer_idx].positions;
                    auto lower = std::lower_bound(positions.begin(), positions.end(), locus.start);
                    auto upper = std::upper_bound(positions.begin(), positions.end(), locus.end);
                    uint32_t count_in_locus = std::distance(lower, upper);
                    if (count_in_locus > max_kmers_in_locus) {
                        max_kmers_in_locus = count_in_locus;
                        signature_kmer_idx = kmer_idx;
                    }
                }
            }
            
            const auto& signature_kmer_data = kmer_db[signature_kmer_idx];

            double entropy = 0.0;
            if (!signature_kmer_data.distances.empty()) {
                std::unordered_map<uint64_t, uint32_t> dist_counts;
                for (uint64_t dist : signature_kmer_data.distances) dist_counts[dist]++;
                for (const auto& pair : dist_counts) {
                    double p = static_cast<double>(pair.second) / signature_kmer_data.distances.size();
                    if (p > 0) entropy -= p * log2(p);
                }
            }

            uint32_t sig_locus_count = 0, sig_max_locus_size = 0;
            if (!signature_kmer_data.positions.empty()) {
                sig_locus_count = 1;
                uint32_t current_locus_size = 1;
                for (size_t p_idx = 1; p_idx < signature_kmer_data.positions.size(); ++p_idx) {
                    if (signature_kmer_data.positions[p_idx] - signature_kmer_data.positions[p_idx - 1] > locus_gap_threshold) {
                        sig_max_locus_size = std::max(sig_max_locus_size, current_locus_size);
                        current_locus_size = 1;
                        sig_locus_count++;
                    } else {
                        current_locus_size++;
                    }
                }
                sig_max_locus_size = std::max(sig_max_locus_size, current_locus_size);
            }

            size_t contig_idx = findContigIndex(locus.start);
            const auto& contig = contig_map[contig_idx];
            uint64_t local_start = locus.start - contig.global_offset;
            uint64_t local_end = locus.end - contig.global_offset + K_SIZE;
            if (local_end > contig.length) local_end = contig.length;

            std::stringstream name_ss;
            name_ss << "locus_" << i+1 << ";"
                    << "kmer=" << uint64ToKmer(signature_kmer_data.kmer_val, K_SIZE) << ";"
                    << "tf=" << signature_kmer_data.tf << ";"
                    << "local_tf=" << max_kmers_in_locus << ";"
                    << "entropy=" << std::fixed << std::setprecision(2) << entropy << ";"
                    << "locus_count=" << sig_locus_count << ";"
                    << "max_locus_size=" << sig_max_locus_size << ";"
                    << "unique_kmers=" << locus.contributing_kmer_indices.size();

            // Filter by minimum locus length
            uint64_t locus_length = local_end - local_start;
            if (locus_length >= min_locus_length) {
                bed_records.push_back({contig.name, local_start, local_end, name_ss.str(), locus_length});
            }
        }
        
        printProgressBar(0.5, "Sorting BED", start_time);
        std::sort(bed_records.begin(), bed_records.end());

        printProgressBar(0.75, "Writing BED file", start_time);
        for (size_t i = 0; i < bed_records.size(); ++i) {
            const auto& rec = bed_records[i];
            out << rec.chrom << "\t" << rec.start << "\t" << rec.end << "\t"
                << rec.name << "\t" << rec.score << "\t" << rec.strand << "\n";
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
    if (argc < 4 || argc > 11) {
        std::cerr << "Usage: " << argv[0] << " <input.bin> <input.fai> <output.bed> [min_tf] [locus_gap] [min_locus_kmers] [inter_overlap] [jaccard_thresh] [min_length] [threads]" << std::endl;
        std::cerr << "  [min_kmer_tf]       - optional. Min. k-mer frequency (default: 100)." << std::endl;
        std::cerr << "  [locus_gap]         - optional. Max. gap within locus (default: 10000)." << std::endl;
        std::cerr << "  [min_locus_kmers]   - optional. Min. number of k-mers in locus (default: 2)." << std::endl;
        std::cerr << "  [inter_overlap]     - optional. Strong overlap threshold for merging (default: 0.9)." << std::endl;
        std::cerr << "  [jaccard_thresh]    - optional. Jaccard index threshold for merging (default: 0.5)." << std::endl;
        std::cerr << "  [min_length]        - optional. Min. locus length for output (default: 100)." << std::endl;
        std::cerr << "  [threads]           - optional. Number of threads (default: all available)." << std::endl;
        return 1;
    }

    std::string input_file = argv[1];
    std::string fai_file = argv[2];
    std::string output_file = argv[3];
    uint32_t min_kmer_tf = 100;
    uint32_t locus_gap_threshold = 10000;
    uint32_t min_locus_kmer_count = 2;
    double intersection_overlap = 0.9;
    double jaccard_threshold = 0.5;
    uint32_t min_locus_length = 100;
    unsigned int num_threads = std::thread::hardware_concurrency();

    if (argc >= 5) { try { min_kmer_tf = std::stoul(argv[4]); } catch (...) { std::cerr << "Error: invalid value for min_kmer_tf." << std::endl; return 1; } }
    if (argc >= 6) { try { locus_gap_threshold = std::stoul(argv[5]); } catch (...) { std::cerr << "Error: invalid value for locus_gap_threshold." << std::endl; return 1; } }
    if (argc >= 7) { try { min_locus_kmer_count = std::stoul(argv[6]); } catch (...) { std::cerr << "Error: invalid value for min_locus_kmers." << std::endl; return 1; } }
    if (argc >= 8) { try { intersection_overlap = std::stod(argv[7]); } catch (...) { std::cerr << "Error: invalid value for intersection_overlap." << std::endl; return 1; } }
    if (argc >= 9) { try { jaccard_threshold = std::stod(argv[8]); } catch (...) { std::cerr << "Error: invalid value for jaccard_threshold." << std::endl; return 1; } }
    if (argc >= 10) { try { min_locus_length = std::stoul(argv[9]); } catch (...) { std::cerr << "Error: invalid value for min_locus_length." << std::endl; return 1; } }
    if (argc >= 11) {
        try {
            num_threads = std::stoul(argv[10]);
            if (num_threads == 0) num_threads = 1;
        } catch (...) { std::cerr << "Error: invalid value for num_threads." << std::endl; return 1; }
    }
    
    try {
        LocusBedGeneratorDetailed merger;
        merger.generate(input_file, fai_file, output_file, min_kmer_tf, locus_gap_threshold, min_locus_kmer_count, intersection_overlap, jaccard_threshold, min_locus_length, num_threads);
        std::cout << "\nAnnotated BED file creation completed successfully." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

