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
#include <thread>
#include <atomic>
#include <mutex>
#include <iomanip>
#include <sstream>
#include <set>
#include <map>

// --- Utilities ---

// Helper function to format time as MM:SS
std::string formatTime(long long seconds) {
    if (seconds < 0) return "--:--";
    long long minutes = seconds / 60;
    seconds %= 60;
    char buffer[32];
    snprintf(buffer, sizeof(buffer), "%02lld:%02lld", minutes, seconds);
    return buffer;
}

// Function to display a progress bar
void printProgressBar(double percentage, const std::string& message, std::chrono::steady_clock::time_point start_time, const std::string& extra_info = "") {
    percentage = std::max(0.0, std::min(1.0, percentage));
    int val = static_cast<int>(percentage * 100);
    int lpad = static_cast<int>(percentage * 40); // Reduced for extra info
    int rpad = 40 - lpad;

    auto now = std::chrono::steady_clock::now();
    double elapsed_sec = std::chrono::duration_cast<std::chrono::duration<double>>(now - start_time).count();
    
    long long eta_sec = -1;
    if (percentage > 0.001) {
        eta_sec = static_cast<long long>(elapsed_sec * (1.0 / percentage - 1.0));
    }
    
    std::string elapsed_str = formatTime(static_cast<long long>(elapsed_sec));
    std::string eta_str = formatTime(eta_sec);

    char buffer[256];
    snprintf(buffer, sizeof(buffer), "\r%3d%% [%.*s%*s] %s | Elapsed: %s | ETA: %s%s%s",
             val, lpad, "========================================", rpad, "",
             message.c_str(), elapsed_str.c_str(), eta_str.c_str(),
             extra_info.empty() ? "" : " | ",
             extra_info.c_str());

    std::cerr << buffer << std::flush;

    if (percentage >= 1.0) {
        std::cerr << std::endl;
    }
}

// Function to convert a 64-bit integer to a k-mer string
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

// Function to convert a k-mer string to a 64-bit integer
uint64_t kmerToUint64(const std::string& kmer_str) {
    uint64_t kmer_val = 0;
    for (char c : kmer_str) {
        kmer_val <<= 2;
        switch (c) {
            case 'A': case 'a': kmer_val |= 0; break;
            case 'C': case 'c': kmer_val |= 1; break;
            case 'G': case 'g': kmer_val |= 2; break;
            case 'T': case 't': kmer_val |= 3; break;
            default: throw std::runtime_error("Invalid nucleotide in k-mer: " + std::string(1, c));
        }
    }
    return kmer_val;
}

// --- Data Structures ---

struct KmerData {
    uint64_t kmer_val;
    uint32_t tf;
    std::vector<uint64_t> positions;
};

struct ContigInfo {
    std::string name;
    uint64_t length;
    uint64_t global_offset;
};

struct KmerContext {
    std::set<uint64_t> upstream_kmers;   // K-mers found upstream
    std::set<uint64_t> downstream_kmers; // K-mers found downstream
};

struct DispersedRepeat {
    uint64_t seed_kmer;
    std::vector<uint64_t> instance_positions;
    KmerContext shared_context;
    double context_similarity_score;
};

struct KmerStats {
    std::string kmer_str;
    uint32_t tf;
    double dist_entropy;
    uint32_t locus_count;
    uint32_t max_locus_size;
};

struct ExclusionRegion {
    std::string chrom;
    uint64_t start;
    uint64_t end;
};

// --- Main Analyzer Class ---

class DispersedRepeatFinder {
public:
    void find(const std::string& input_binary, const std::string& input_cidx, const std::string& output_bed,
              const std::string& distance_analysis_file, const std::string& exclusion_bed_file,
              uint32_t min_tf, uint32_t window_size, double min_context_similarity,
              uint32_t min_instances, double max_entropy, unsigned int num_threads) {
        
        // Load contig information
        loadContigIndex(input_cidx);
        
        // Load k-mer statistics if provided
        std::unordered_map<std::string, KmerStats> kmer_stats;
        if (!distance_analysis_file.empty()) {
            loadKmerStats(distance_analysis_file, kmer_stats);
        }
        
        // Load exclusion regions if provided
        std::vector<ExclusionRegion> exclusion_regions;
        if (!exclusion_bed_file.empty()) {
            loadExclusionRegions(exclusion_bed_file, exclusion_regions);
        }
        
        // Load k-mer data
        std::vector<KmerData> kmer_data;
        loadKmerData(input_binary, kmer_data, min_tf);
        
        if (kmer_data.empty()) {
            std::cerr << "No k-mers found with frequency >= " << min_tf << std::endl;
            return;
        }
        
        // Filter k-mers based on entropy if stats are available
        if (!kmer_stats.empty()) {
            filterByEntropy(kmer_data, kmer_stats, max_entropy);
        }
        
        // Build position index for fast context lookup
        buildPositionIndex(kmer_data);
        
        // Find dispersed repeats
        std::vector<DispersedRepeat> dispersed_repeats;
        findDispersedRepeats(kmer_data, dispersed_repeats, exclusion_regions, window_size, 
                           min_context_similarity, min_instances, num_threads);
        
        // Write results
        writeBedFile(output_bed, dispersed_repeats);
    }

private:
    std::vector<ContigInfo> contig_map;
    std::unordered_map<uint64_t, std::vector<uint64_t>> position_to_kmers; // position -> list of k-mers at that position
    const int K_SIZE = 13;
    
    void loadContigIndex(const std::string& filename) {
        auto start_time = std::chrono::steady_clock::now();
        std::ifstream in(filename);
        if (!in) throw std::runtime_error("Cannot open contig index file: " + filename);
        
        std::string line;
        while (std::getline(in, line)) {
            std::istringstream iss(line);
            ContigInfo info;
            iss >> info.name >> info.length >> info.global_offset;
            contig_map.push_back(info);
        }
        
        printProgressBar(1.0, "Loading contig index", start_time);
    }
    
    void loadKmerStats(const std::string& filename, std::unordered_map<std::string, KmerStats>& kmer_stats) {
        auto start_time = std::chrono::steady_clock::now();
        std::ifstream in(filename);
        if (!in) {
            std::cerr << "Warning: Cannot open distance analysis file: " << filename << std::endl;
            return;
        }
        
        std::string line;
        // Skip header
        std::getline(in, line);
        
        while (std::getline(in, line)) {
            std::istringstream iss(line);
            KmerStats stats;
            std::string distances_summary;
            
            if (iss >> stats.kmer_str >> stats.tf >> stats.dist_entropy 
                >> stats.locus_count >> stats.max_locus_size) {
                kmer_stats[stats.kmer_str] = stats;
            }
        }
        
        printProgressBar(1.0, "Loading k-mer statistics", start_time);
        std::cerr << "Loaded statistics for " << kmer_stats.size() << " k-mers" << std::endl;
    }
    
    void loadExclusionRegions(const std::string& filename, std::vector<ExclusionRegion>& exclusion_regions) {
        auto start_time = std::chrono::steady_clock::now();
        std::ifstream in(filename);
        if (!in) {
            std::cerr << "Warning: Cannot open exclusion BED file: " << filename << std::endl;
            return;
        }
        
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#' || line.substr(0, 5) == "track") continue;
            
            std::istringstream iss(line);
            ExclusionRegion region;
            if (iss >> region.chrom >> region.start >> region.end) {
                exclusion_regions.push_back(region);
            }
        }
        
        printProgressBar(1.0, "Loading exclusion regions", start_time);
        std::cerr << "Loaded " << exclusion_regions.size() << " exclusion regions" << std::endl;
    }
    
    void filterByEntropy(std::vector<KmerData>& kmer_data, 
                        const std::unordered_map<std::string, KmerStats>& kmer_stats,
                        double max_entropy) {
        auto start_time = std::chrono::steady_clock::now();
        std::vector<KmerData> filtered;
        
        for (const auto& kmer : kmer_data) {
            std::string kmer_str = uint64ToKmer(kmer.kmer_val, K_SIZE);
            auto it = kmer_stats.find(kmer_str);
            
            if (it != kmer_stats.end()) {
                // High entropy = more dispersed, low entropy = more clustered
                if (it->second.dist_entropy >= max_entropy) {
                    filtered.push_back(kmer);
                }
            } else {
                // Keep k-mers without stats
                filtered.push_back(kmer);
            }
        }
        
        printProgressBar(1.0, "Filtering by entropy", start_time);
        std::cerr << "Filtered from " << kmer_data.size() << " to " << filtered.size() 
                  << " k-mers (entropy >= " << max_entropy << ")" << std::endl;
        
        kmer_data = std::move(filtered);
    }
    
    void loadKmerData(const std::string& filename, std::vector<KmerData>& kmer_data, uint32_t min_tf) {
        auto start_time = std::chrono::steady_clock::now();
        std::ifstream in(filename, std::ios::binary);
        if (!in) throw std::runtime_error("Cannot open binary file: " + filename);
        
        // Check magic number and version
        uint32_t magic;
        in.read(reinterpret_cast<char*>(&magic), sizeof(magic));
        if (magic != 0x524D4B41) { // "AKMR"
            // Old format without header
            in.seekg(0);
        } else {
            uint32_t version;
            in.read(reinterpret_cast<char*>(&version), sizeof(version));
            if (version != 1) {
                throw std::runtime_error("Unsupported file version");
            }
        }
        
        uint64_t total_kmers;
        in.read(reinterpret_cast<char*>(&total_kmers), sizeof(total_kmers));
        
        printProgressBar(0.0, "Loading k-mer data", start_time);
        
        for (uint64_t i = 0; i < total_kmers; ++i) {
            KmerData kmer;
            in.read(reinterpret_cast<char*>(&kmer.kmer_val), sizeof(kmer.kmer_val));
            in.read(reinterpret_cast<char*>(&kmer.tf), sizeof(kmer.tf));
            
            if (kmer.tf < min_tf) break; // Data is sorted by frequency
            
            uint64_t pos_size, dist_size;
            in.read(reinterpret_cast<char*>(&pos_size), sizeof(pos_size));
            kmer.positions.resize(pos_size);
            in.read(reinterpret_cast<char*>(kmer.positions.data()), pos_size * sizeof(uint64_t));
            
            // Skip distances
            in.read(reinterpret_cast<char*>(&dist_size), sizeof(dist_size));
            in.seekg(dist_size * sizeof(uint64_t), std::ios::cur);
            
            kmer_data.push_back(kmer);
            
            if (i % 1000 == 0) {
                printProgressBar(static_cast<double>(i) / total_kmers, "Loading k-mer data", start_time);
            }
        }
        
        printProgressBar(1.0, "Loading k-mer data", start_time);
    }
    
    void buildPositionIndex(const std::vector<KmerData>& kmer_data) {
        auto start_time = std::chrono::steady_clock::now();
        printProgressBar(0.0, "Building position index", start_time);
        
        for (size_t i = 0; i < kmer_data.size(); ++i) {
            const auto& kmer = kmer_data[i];
            for (uint64_t pos : kmer.positions) {
                position_to_kmers[pos].push_back(kmer.kmer_val);
            }
            
            if (i % 1000 == 0) {
                printProgressBar(static_cast<double>(i) / kmer_data.size(), "Building position index", start_time);
            }
        }
        
        printProgressBar(1.0, "Building position index", start_time);
    }
    
    KmerContext getKmerContext(uint64_t position, uint32_t window_size) {
        KmerContext context;
        
        // Check upstream
        uint64_t upstream_start = (position > window_size) ? position - window_size : 0;
        for (uint64_t pos = upstream_start; pos < position; ++pos) {
            auto it = position_to_kmers.find(pos);
            if (it != position_to_kmers.end()) {
                for (uint64_t kmer : it->second) {
                    context.upstream_kmers.insert(kmer);
                }
            }
        }
        
        // Check downstream
        uint64_t downstream_end = position + K_SIZE + window_size;
        for (uint64_t pos = position + K_SIZE; pos < downstream_end; ++pos) {
            auto it = position_to_kmers.find(pos);
            if (it != position_to_kmers.end()) {
                for (uint64_t kmer : it->second) {
                    context.downstream_kmers.insert(kmer);
                }
            }
        }
        
        return context;
    }
    
    double calculateContextSimilarity(const KmerContext& ctx1, const KmerContext& ctx2) {
        if (ctx1.upstream_kmers.empty() && ctx1.downstream_kmers.empty() &&
            ctx2.upstream_kmers.empty() && ctx2.downstream_kmers.empty()) {
            return 0.0;
        }
        
        // Calculate Jaccard similarity for upstream
        std::set<uint64_t> upstream_intersection;
        std::set_intersection(ctx1.upstream_kmers.begin(), ctx1.upstream_kmers.end(),
                            ctx2.upstream_kmers.begin(), ctx2.upstream_kmers.end(),
                            std::inserter(upstream_intersection, upstream_intersection.begin()));
        
        std::set<uint64_t> upstream_union;
        std::set_union(ctx1.upstream_kmers.begin(), ctx1.upstream_kmers.end(),
                      ctx2.upstream_kmers.begin(), ctx2.upstream_kmers.end(),
                      std::inserter(upstream_union, upstream_union.begin()));
        
        double upstream_jaccard = upstream_union.empty() ? 0.0 : 
            static_cast<double>(upstream_intersection.size()) / upstream_union.size();
        
        // Calculate Jaccard similarity for downstream
        std::set<uint64_t> downstream_intersection;
        std::set_intersection(ctx1.downstream_kmers.begin(), ctx1.downstream_kmers.end(),
                            ctx2.downstream_kmers.begin(), ctx2.downstream_kmers.end(),
                            std::inserter(downstream_intersection, downstream_intersection.begin()));
        
        std::set<uint64_t> downstream_union;
        std::set_union(ctx1.downstream_kmers.begin(), ctx1.downstream_kmers.end(),
                      ctx2.downstream_kmers.begin(), ctx2.downstream_kmers.end(),
                      std::inserter(downstream_union, downstream_union.begin()));
        
        double downstream_jaccard = downstream_union.empty() ? 0.0 :
            static_cast<double>(downstream_intersection.size()) / downstream_union.size();
        
        // Average of upstream and downstream similarities
        return (upstream_jaccard + downstream_jaccard) / 2.0;
    }
    
    bool isInExclusionRegion(uint64_t global_pos, const std::vector<ExclusionRegion>& exclusion_regions) {
        if (exclusion_regions.empty()) return false;
        
        size_t contig_idx = findContigIndex(global_pos);
        const auto& contig = contig_map[contig_idx];
        uint64_t local_pos = global_pos - contig.global_offset;
        
        for (const auto& region : exclusion_regions) {
            if (region.chrom == contig.name && local_pos >= region.start && local_pos < region.end) {
                return true;
            }
        }
        return false;
    }
    
    void findDispersedRepeats(const std::vector<KmerData>& kmer_data,
                             std::vector<DispersedRepeat>& dispersed_repeats,
                             const std::vector<ExclusionRegion>& exclusion_regions,
                             uint32_t window_size, double min_context_similarity,
                             uint32_t min_instances, unsigned int num_threads) {
        auto start_time = std::chrono::steady_clock::now();
        std::mutex result_mutex;
        
        // Process k-mers in parallel
        std::atomic<size_t> next_kmer_idx{0};
        std::vector<std::thread> threads;
        
        auto worker = [&]() {
            while (true) {
                size_t idx = next_kmer_idx.fetch_add(1);
                if (idx >= kmer_data.size()) break;
                
                const auto& kmer = kmer_data[idx];
                if (kmer.positions.size() < min_instances) continue;
                
                // Get contexts for all instances, excluding those in tandem repeat regions
                std::vector<KmerContext> contexts;
                std::vector<uint64_t> valid_positions;
                for (uint64_t pos : kmer.positions) {
                    if (!isInExclusionRegion(pos, exclusion_regions)) {
                        contexts.push_back(getKmerContext(pos, window_size));
                        valid_positions.push_back(pos);
                    }
                }
                
                if (valid_positions.size() < min_instances) continue;
                
                // Find groups of instances with similar contexts
                std::vector<bool> used(contexts.size(), false);
                
                for (size_t i = 0; i < contexts.size(); ++i) {
                    if (used[i]) continue;
                    
                    DispersedRepeat repeat;
                    repeat.seed_kmer = kmer.kmer_val;
                    repeat.instance_positions.push_back(valid_positions[i]);
                    repeat.shared_context = contexts[i];
                    
                    // Find other instances with similar context
                    for (size_t j = i + 1; j < contexts.size(); ++j) {
                        if (used[j]) continue;
                        
                        double similarity = calculateContextSimilarity(contexts[i], contexts[j]);
                        if (similarity >= min_context_similarity) {
                            repeat.instance_positions.push_back(valid_positions[j]);
                            used[j] = true;
                            
                            // Update shared context (intersection)
                            std::set<uint64_t> new_upstream;
                            std::set_intersection(repeat.shared_context.upstream_kmers.begin(),
                                                repeat.shared_context.upstream_kmers.end(),
                                                contexts[j].upstream_kmers.begin(),
                                                contexts[j].upstream_kmers.end(),
                                                std::inserter(new_upstream, new_upstream.begin()));
                            repeat.shared_context.upstream_kmers = new_upstream;
                            
                            std::set<uint64_t> new_downstream;
                            std::set_intersection(repeat.shared_context.downstream_kmers.begin(),
                                                repeat.shared_context.downstream_kmers.end(),
                                                contexts[j].downstream_kmers.begin(),
                                                contexts[j].downstream_kmers.end(),
                                                std::inserter(new_downstream, new_downstream.begin()));
                            repeat.shared_context.downstream_kmers = new_downstream;
                        }
                    }
                    
                    if (repeat.instance_positions.size() >= min_instances) {
                        // Calculate average similarity score
                        double total_similarity = 0.0;
                        int count = 0;
                        for (size_t j = 0; j < repeat.instance_positions.size(); ++j) {
                            for (size_t k = j + 1; k < repeat.instance_positions.size(); ++k) {
                                total_similarity += calculateContextSimilarity(
                                    getKmerContext(repeat.instance_positions[j], window_size),
                                    getKmerContext(repeat.instance_positions[k], window_size)
                                );
                                count++;
                            }
                        }
                        repeat.context_similarity_score = count > 0 ? total_similarity / count : 0.0;
                        
                        std::lock_guard<std::mutex> lock(result_mutex);
                        dispersed_repeats.push_back(repeat);
                    }
                    
                    used[i] = true;
                }
            }
        };
        
        for (unsigned int i = 0; i < num_threads; ++i) {
            threads.emplace_back(worker);
        }
        
        // Monitor progress
        size_t total = kmer_data.size();
        while (next_kmer_idx < total) {
            double progress = static_cast<double>(next_kmer_idx) / total;
            printProgressBar(progress, "Finding dispersed repeats", start_time);
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        
        for (auto& t : threads) {
            t.join();
        }
        
        printProgressBar(1.0, "Finding dispersed repeats", start_time);
    }
    
    size_t findContigIndex(uint64_t global_pos) {
        auto it = std::upper_bound(contig_map.begin(), contig_map.end(), global_pos,
            [](uint64_t pos, const ContigInfo& info) { return pos < info.global_offset; });
        if (it == contig_map.begin()) return 0;
        return std::distance(contig_map.begin(), it) - 1;
    }
    
    void writeBedFile(const std::string& filename, const std::vector<DispersedRepeat>& dispersed_repeats) {
        auto start_time = std::chrono::steady_clock::now();
        std::ofstream out(filename);
        if (!out) throw std::runtime_error("Cannot open output file: " + filename);
        
        out << "track name=\"DispersedRepeats\" description=\"Dispersed repeats with similar k-mer contexts\"\n";
        
        printProgressBar(0.0, "Writing BED file", start_time);
        
        for (size_t i = 0; i < dispersed_repeats.size(); ++i) {
            const auto& repeat = dispersed_repeats[i];
            std::string seed_kmer_str = uint64ToKmer(repeat.seed_kmer, K_SIZE);
            
            // Write each instance as a separate BED entry
            for (size_t j = 0; j < repeat.instance_positions.size(); ++j) {
                uint64_t pos = repeat.instance_positions[j];
                size_t contig_idx = findContigIndex(pos);
                const auto& contig = contig_map[contig_idx];
                uint64_t local_pos = pos - contig.global_offset;
                
                // Create info string with context k-mers
                std::stringstream info;
                info << "repeat_" << i+1 << "_instance_" << j+1 << ";"
                     << "seed=" << seed_kmer_str << ";"
                     << "instances=" << repeat.instance_positions.size() << ";"
                     << "similarity=" << std::fixed << std::setprecision(3) << repeat.context_similarity_score << ";"
                     << "upstream_context=" << repeat.shared_context.upstream_kmers.size() << ";"
                     << "downstream_context=" << repeat.shared_context.downstream_kmers.size();
                
                out << contig.name << "\t"
                    << local_pos << "\t"
                    << local_pos + K_SIZE << "\t"
                    << info.str() << "\t"
                    << static_cast<int>(repeat.context_similarity_score * 1000) << "\t"
                    << ".\n";
            }
            
            if (i % 100 == 0) {
                printProgressBar(static_cast<double>(i) / dispersed_repeats.size(), "Writing BED file", start_time);
            }
        }
        
        printProgressBar(1.0, "Writing BED file", start_time);
    }
};

int main(int argc, char* argv[]) {
    if (argc < 4 || argc > 12) {
        std::cerr << "Usage: " << argv[0] << " <input.bin> <contigs.cidx> <output.bed> [options]\n";
        std::cerr << "Required arguments:\n";
        std::cerr << "  <input.bin>      - Binary file from kmer_analyzer\n";
        std::cerr << "  <contigs.cidx>   - Contig index file from kmer_analyzer\n";
        std::cerr << "  <output.bed>     - Output BED file with dispersed repeats\n";
        std::cerr << "\nOptional arguments:\n";
        std::cerr << "  [distance.tsv]   - Distance analysis file for entropy filtering (use 'none' to skip)\n";
        std::cerr << "  [exclusion.bed]  - BED file with regions to exclude (use 'none' to skip)\n";
        std::cerr << "  [min_tf]         - Minimum k-mer frequency (default: 100)\n";
        std::cerr << "  [window]         - Context window size in bp (default: 500)\n";
        std::cerr << "  [similarity]     - Minimum context similarity (0-1, default: 0.7)\n";
        std::cerr << "  [min_inst]       - Minimum instances per repeat (default: 5)\n";
        std::cerr << "  [max_entropy]    - Maximum entropy for seed selection (default: 1.5)\n";
        std::cerr << "  [threads]        - Number of threads (default: all available)\n";
        return 1;
    }
    
    std::string input_file = argv[1];
    std::string cidx_file = argv[2];
    std::string output_file = argv[3];
    
    // Optional file arguments
    std::string distance_file = "";
    std::string exclusion_file = "";
    
    // Default values
    uint32_t min_tf = 100;
    uint32_t window_size = 500;
    double min_similarity = 0.7;
    uint32_t min_instances = 5;
    double max_entropy = 1.5;
    unsigned int num_threads = std::thread::hardware_concurrency();
    
    int arg_idx = 4;
    
    // Parse optional files
    if (argc > arg_idx && std::string(argv[arg_idx]) != "none") {
        distance_file = argv[arg_idx];
    }
    arg_idx++;
    
    if (argc > arg_idx && std::string(argv[arg_idx]) != "none") {
        exclusion_file = argv[arg_idx];
    }
    arg_idx++;
    
    // Parse numeric parameters
    if (argc > arg_idx) {
        try { min_tf = std::stoul(argv[arg_idx]); }
        catch (...) { std::cerr << "Error: invalid min_tf\n"; return 1; }
        arg_idx++;
    }
    if (argc > arg_idx) {
        try { window_size = std::stoul(argv[arg_idx]); }
        catch (...) { std::cerr << "Error: invalid window size\n"; return 1; }
        arg_idx++;
    }
    if (argc > arg_idx) {
        try { min_similarity = std::stod(argv[arg_idx]); }
        catch (...) { std::cerr << "Error: invalid similarity threshold\n"; return 1; }
        arg_idx++;
    }
    if (argc > arg_idx) {
        try { min_instances = std::stoul(argv[arg_idx]); }
        catch (...) { std::cerr << "Error: invalid min_instances\n"; return 1; }
        arg_idx++;
    }
    if (argc > arg_idx) {
        try { max_entropy = std::stod(argv[arg_idx]); }
        catch (...) { std::cerr << "Error: invalid max_entropy\n"; return 1; }
        arg_idx++;
    }
    if (argc > arg_idx) {
        try { num_threads = std::stoul(argv[arg_idx]); }
        catch (...) { std::cerr << "Error: invalid thread count\n"; return 1; }
        arg_idx++;
    }
    
    try {
        DispersedRepeatFinder finder;
        finder.find(input_file, cidx_file, output_file, distance_file, exclusion_file,
                   min_tf, window_size, min_similarity, min_instances, max_entropy, num_threads);
        std::cout << "\nDispersed repeat finding completed successfully.\n";
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << "\n";
        return 1;
    }
    
    return 0;
}