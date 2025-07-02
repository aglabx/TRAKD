#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <chrono>
#include <iomanip> // для std::setprecision
#include <sstream> // для std::stringstream
#include <thread>
#include <atomic>
#include <mutex>
#include <cmath> // для log2

// --- Утилиты ---

// Вспомогательная функция для форматирования времени в MM:SS
std::string formatTime(long long seconds) {
    if (seconds < 0) return "--:--";
    long long minutes = seconds / 60;
    seconds %= 60;
    char buffer[32];
    snprintf(buffer, sizeof(buffer), "%02lld:%02lld", minutes, seconds);
    return buffer;
}

// Функция для отображения прогресс-бара
void printProgressBar(double percentage, const std::string& message, std::chrono::steady_clock::time_point start_time, const std::string& extra_info = "") {
    percentage = std::max(0.0, std::min(1.0, percentage));
    int val = static_cast<int>(percentage * 100);
    int lpad = static_cast<int>(percentage * 40); // Уменьшено для доп. информации
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

// Функция для преобразования 64-битного числа в k-mer строку
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

// --- Основной класс анализатора ---

struct KmerRawData {
    uint64_t kmer_val;
    uint32_t kmer_tf;
    std::vector<uint32_t> positions;
    std::vector<uint32_t> distances;
};

struct KmerAnalysisResult {
    uint64_t original_index;
    std::string kmer_str;
    uint32_t kmer_tf;
    double dist_entropy;
    uint32_t locus_count;
    uint32_t max_locus_size;
    std::string summary;
};

class DistanceAnalyzer {
public:
    void analyze(const std::string& input_binary_full, const std::string& output_tsv, int k_size, uint32_t min_kmer_tf, double min_dist_fraction, uint32_t locus_gap_threshold, unsigned int num_threads) {
        
        // --- Этап 1: Чтение всех данных в память ---
        std::vector<KmerRawData> jobs;
        readAllJobs(input_binary_full, jobs, min_kmer_tf);
        
        if (jobs.empty()) {
            std::cerr << "Нет данных для анализа после применения фильтров." << std::endl;
            // Создаем пустой файл с заголовком
            std::ofstream out(output_tsv);
            out << "kmer\tkmer_tf\tdist_entropy\tlocus_count\tmax_locus_size\tdistances_summary\n";
            return;
        }

        // --- Этап 2: Параллельный анализ ---
        std::vector<KmerAnalysisResult> results;
        processJobsParallel(jobs, results, k_size, min_dist_fraction, locus_gap_threshold, num_threads);

        // --- Этап 3: Запись результатов ---
        writeResultsToFile(output_tsv, results);
    }

private:
    void readAllJobs(const std::string& filename, std::vector<KmerRawData>& jobs, uint32_t min_kmer_tf) {
        auto start_time = std::chrono::steady_clock::now();
        std::ifstream in(filename, std::ios::binary);
        if (!in) throw std::runtime_error("Невозможно открыть полный бинарный файл для чтения: " + filename);
        
        uint64_t total_kmers_in_file;
        in.read(reinterpret_cast<char*>(&total_kmers_in_file), sizeof(total_kmers_in_file));
        if (!in || total_kmers_in_file == 0) {
            printProgressBar(1.0, "Чтение заданий", start_time);
            return;
        }
        
        printProgressBar(0.0, "Чтение заданий", start_time);
        
        for (uint64_t i = 0; i < total_kmers_in_file; ++i) {
            KmerRawData job;
            in.read(reinterpret_cast<char*>(&job.kmer_val), sizeof(job.kmer_val));
            in.read(reinterpret_cast<char*>(&job.kmer_tf), sizeof(job.kmer_tf));

            if (job.kmer_tf < min_kmer_tf) break;

            uint64_t pos_size, dist_size;
            in.read(reinterpret_cast<char*>(&pos_size), sizeof(pos_size));
            job.positions.resize(pos_size);
            in.read(reinterpret_cast<char*>(job.positions.data()), pos_size * sizeof(uint32_t));

            in.read(reinterpret_cast<char*>(&dist_size), sizeof(dist_size));
            job.distances.resize(dist_size);
            in.read(reinterpret_cast<char*>(job.distances.data()), dist_size * sizeof(uint32_t));

            if (!in) throw std::runtime_error("Ошибка при чтении бинарного файла.");
            jobs.push_back(std::move(job));
            
            if (i % 1000 == 0) {
                printProgressBar(static_cast<double>(i + 1) / total_kmers_in_file, "Чтение заданий", start_time);
            }
        }
        printProgressBar(1.0, "Чтение заданий", start_time);
    }

    void worker_thread_function(
        const std::vector<KmerRawData>& jobs,
        std::vector<KmerAnalysisResult>& results,
        std::atomic<size_t>& next_job_idx,
        int k_size,
        double min_dist_fraction,
        uint32_t locus_gap_threshold
    ) {
        while (true) {
            size_t i = next_job_idx.fetch_add(1);
            if (i >= jobs.size()) break;

            const auto& job = jobs[i];
            
            // --- Анализ дистанций ---
            std::unordered_map<uint32_t, uint32_t> local_dist_counts;
            for (uint32_t dist : job.distances) {
                local_dist_counts[dist]++;
            }
            
            double entropy = 0.0;
            if (!job.distances.empty()) {
                for (const auto& pair : local_dist_counts) {
                    double p = static_cast<double>(pair.second) / job.distances.size();
                    if (p > 0) { // log2(0) не определен
                        entropy -= p * log2(p);
                    }
                }
            }

            std::vector<std::pair<uint32_t, uint32_t>> sorted_local_dists;
            for (const auto& pair : local_dist_counts) sorted_local_dists.push_back(pair);
            std::sort(sorted_local_dists.begin(), sorted_local_dists.end(), 
                [](const auto& a, const auto& b) { return a.second > b.second; });

            std::stringstream summary_ss;
            int distances_to_write = std::min(100, static_cast<int>(sorted_local_dists.size()));
            bool first_item = true;
            for (int j = 0; j < distances_to_write; ++j) {
                const auto& pair = sorted_local_dists[j];
                double fraction = (job.distances.size() > 0) ? static_cast<double>(pair.second) / job.distances.size() : 0.0;
                if (fraction >= min_dist_fraction) {
                    if (!first_item) summary_ss << "|";
                    summary_ss << pair.first << " " << pair.second << " " << std::fixed << std::setprecision(6) << fraction;
                    first_item = false;
                }
            }

            // --- Анализ локусов ---
            uint32_t locus_count = 0;
            uint32_t max_locus_size = 0;
            if (!job.positions.empty()) {
                locus_count = 1;
                uint32_t current_locus_size = 1;
                for (size_t p_idx = 1; p_idx < job.positions.size(); ++p_idx) {
                    if (job.positions[p_idx] - job.positions[p_idx - 1] > locus_gap_threshold) {
                        max_locus_size = std::max(max_locus_size, current_locus_size);
                        current_locus_size = 1;
                        locus_count++;
                    } else {
                        current_locus_size++;
                    }
                }
                max_locus_size = std::max(max_locus_size, current_locus_size);
            }

            results[i] = {i, uint64ToKmer(job.kmer_val, k_size), job.kmer_tf, entropy, locus_count, max_locus_size, summary_ss.str()};
        }
    }

    void processJobsParallel(
        const std::vector<KmerRawData>& jobs,
        std::vector<KmerAnalysisResult>& results,
        int k_size,
        double min_dist_fraction,
        uint32_t locus_gap_threshold,
        unsigned int num_threads
    ) {
        auto start_time = std::chrono::steady_clock::now();
        std::string msg = "Анализ (" + std::to_string(num_threads) + " п.)";
        results.resize(jobs.size());

        std::vector<std::thread> threads;
        std::atomic<size_t> next_job_idx{0};
        
        for (unsigned int i = 0; i < num_threads; ++i) {
            threads.emplace_back(&DistanceAnalyzer::worker_thread_function, this, std::ref(jobs), std::ref(results), std::ref(next_job_idx), k_size, min_dist_fraction, locus_gap_threshold);
        }
        
        while(true) {
            size_t processed_count = next_job_idx.load();
            if (processed_count > jobs.size()) processed_count = jobs.size();
            std::string extra_info = "[" + std::to_string(processed_count) + "/" + std::to_string(jobs.size()) + "]";
            printProgressBar(static_cast<double>(processed_count) / jobs.size(), msg, start_time, extra_info);
            if (processed_count >= jobs.size()) break;
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }

        for (auto& t : threads) {
            t.join();
        }
    }

    void writeResultsToFile(const std::string& filename, std::vector<KmerAnalysisResult>& results) {
        auto start_time = std::chrono::steady_clock::now();
        std::ofstream out(filename);
        if (!out) throw std::runtime_error("Невозможно открыть TSV файл для записи: " + filename);

        out << "kmer\tkmer_tf\tdist_entropy\tlocus_count\tmax_locus_size\tdistances_summary\n";
        
        if (results.empty()) {
            printProgressBar(1.0, "Запись результатов", start_time);
            return;
        }

        printProgressBar(0.0, "Запись результатов", start_time);
        for(size_t i = 0; i < results.size(); ++i) {
            const auto& res = results[i];
            out << res.kmer_str << "\t" 
                << res.kmer_tf << "\t" 
                << std::fixed << std::setprecision(4) << res.dist_entropy << "\t"
                << res.locus_count << "\t"
                << res.max_locus_size << "\t"
                << res.summary << "\n";
            if (i % 1000 == 0 || i == results.size() - 1) {
                printProgressBar(static_cast<double>(i + 1) / results.size(), "Запись результатов", start_time);
            }
        }
        printProgressBar(1.0, "Запись результатов", start_time);
    }
};

int main(int argc, char* argv[]) {
    if (argc < 3 || argc > 7) {
        std::cerr << "Использование: " << argv[0] << " <input.bin> <output.tsv> [min_tf] [min_frac] [locus_gap] [threads]" << std::endl;
        std::cerr << "  <input_full.bin>       - бинарный файл от kmer_analyzer." << std::endl;
        std::cerr << "  <output_analysis.tsv>  - TSV файл с детальным анализом." << std::endl;
        std::cerr << "  [min_kmer_tf]          - необязательный. Мин. частота k-мера (по умолч.: 0)." << std::endl;
        std::cerr << "  [min_dist_fraction]    - необязательный. Мин. доля дистанции (по умолч.: 0.0)." << std::endl;
        std::cerr << "  [locus_gap_threshold]  - необязательный. Макс. разрыв внутри локуса (по умолч.: 10000)." << std::endl;
        std::cerr << "  [num_threads]          - необязательный. Количество потоков (по умолч.: все доступные)." << std::endl;
        return 1;
    }

    std::string input_file = argv[1];
    std::string output_file = argv[2];
    uint32_t min_kmer_tf = 0;
    double min_dist_fraction = 0.0;
    uint32_t locus_gap_threshold = 10000;
    unsigned int num_threads = std::thread::hardware_concurrency();

    if (argc >= 4) { try { min_kmer_tf = std::stoul(argv[3]); } catch (...) { std::cerr << "Ошибка: неверное значение для min_kmer_tf." << std::endl; return 1; } }
    if (argc >= 5) { try { min_dist_fraction = std::stod(argv[4]); } catch (...) { std::cerr << "Ошибка: неверное значение для min_dist_fraction." << std::endl; return 1; } }
    if (argc >= 6) { try { locus_gap_threshold = std::stoul(argv[5]); } catch (...) { std::cerr << "Ошибка: неверное значение для locus_gap_threshold." << std::endl; return 1; } }
    if (argc >= 7) {
        try {
            num_threads = std::stoul(argv[6]);
            if (num_threads == 0) { num_threads = 1; }
        } catch (...) { std::cerr << "Ошибка: неверное значение для num_threads." << std::endl; return 1; }
    }

    const int K_SIZE = 13;
    
    try {
        DistanceAnalyzer analyzer;
        analyzer.analyze(input_file, output_file, K_SIZE, min_kmer_tf, min_dist_fraction, locus_gap_threshold, num_threads);
        std::cout << "\nДетальный анализ дистанций успешно завершен." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "\nОшибка: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

