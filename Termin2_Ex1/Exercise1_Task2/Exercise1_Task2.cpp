
#include <iostream>
#include <vector>
#include <random>
#include <array>
#include <limits>
#include <cmath>
#include <chrono>
#include <omp.h>
#include <fstream>

constexpr size_t n_samples = 1000000;
constexpr size_t n_dimensions = 3;
constexpr size_t n_centroids = 10;
constexpr size_t max_iterations = 10000;
const std::pair<float, float> space = { 0.0f, 100.0f };

 


using Point = std::array<float, n_dimensions>;
using Samples = std::vector<Point>;
using Centroids = std::vector<Point>;
using CentroidAssignment = std::vector<int>;

void generate_initial_state(Samples& samples, Centroids& centroids, CentroidAssignment& centroid_assignment) {
    std::mt19937 rng(42);
    std::uniform_real_distribution<float> dist(space.first, space.second);
    for (int i = 0; i < n_samples; ++i) {
        for (int d = 0; d < n_dimensions; ++d)
            samples[i][d] = dist(rng);
    }

    for (int i = 0; i < n_centroids; ++i) {
        for (int d = 0; d < n_dimensions; ++d)
            centroids[i][d] = dist(rng);
	}

	centroid_assignment.assign(n_samples, -1);
}

void assign_samples_to_centroid(const Samples& samples, const Centroids& centroids, CentroidAssignment& centroid_assignment, int num_threads) {
    
    #pragma omp parallel for num_threads(num_threads) 
    for (int i = 0; i < n_samples; ++i) {
        float min_distance = std::numeric_limits<float>::max();
        int closest_centroid = -1;
        for (int j = 0; j < n_centroids; ++j) {
            float distance = 0.0f;
            for (int d = 0; d < n_dimensions; ++d) {
                float diff = samples[i][d] - centroids[j][d];
                distance += diff * diff;
            }
            //distance = std::sqrt(distance);
            if (distance < min_distance) {
                min_distance = distance;
                closest_centroid = j;
            }
        }
        centroid_assignment[i] = closest_centroid;
    }
    
}

Centroids update_centroids(const Samples& samples, Centroids& centroids, const CentroidAssignment& centroid_assignment) {
    
    Centroids new_centroids(n_centroids);
    for (auto& p : new_centroids) p.fill(0.0f);
    std::array<int, n_centroids> counts = {};

    // sum all assigned samples for each centroid
    for (size_t i = 0; i < n_centroids; ++i) {
        for (size_t j = 0; j < n_samples; ++j) {
            if (centroid_assignment[j] == static_cast<int>(i)) {

                for (size_t d = 0; d < n_dimensions; ++d)
                    new_centroids[i][d] += samples[j][d];
                counts[i]++;
            }
		}

    }

    // assign new centroid positions
    for (size_t j = 0; j < n_centroids; ++j) {
        if (counts[j] > 0) {
            for (size_t d = 0; d < n_dimensions; ++d)
                new_centroids[j][d] /= static_cast<float>(counts[j]);  // divide by number of samples for mean value
        }
        else {
			new_centroids[j] = centroids[j];
        }
        // else: no samples assigned -> Centroid unchanged
    }

    return new_centroids;
}

int calculate_cluster(Samples& samples, Centroids& centroids, CentroidAssignment& centroid_assignment, int num_threads) {
	std::cout << "Calculating Clusters using " << num_threads << " Threads." << std::endl;
    for (size_t i = 0; i <= max_iterations; ++i) {
        std::cout << "Start Iteration " << i << std::endl;
        assign_samples_to_centroid(samples, centroids, centroid_assignment, num_threads);
        Centroids new_centroids = update_centroids(samples, centroids, centroid_assignment);

        if (new_centroids == centroids) {
            std::cout << "Converged after " << i << " iterations." << std::endl;
            return i;
        }

    }
}


int main()
{
#ifdef _OPENMP
    std::cout << "_OPENMP defined, value: " << _OPENMP << '\n';
#else
    std::cout << "_OPENMP NOT defined — enable __OpenMP Support__ (/openmp) in Project Properties\n";
#endif

    Samples samples(n_samples);
	Centroids centroids(n_centroids);
	CentroidAssignment centroid_assignment(n_samples);

    generate_initial_state(samples, centroids, centroid_assignment);

    Samples samples_serial = samples;
	Centroids centroids_serial = centroids;
	CentroidAssignment centroid_assignment_serial = centroid_assignment;

    Samples samples_parallel_2 = samples;
    Centroids centroids_parallel_2 = centroids;
	CentroidAssignment centroid_assignment_parallel_2 = centroid_assignment;

    Samples samples_parallel_4 = samples;
    Centroids centroids_parallel_4 = centroids;
    CentroidAssignment centroid_assignment_parallel_4 = centroid_assignment;

    Samples samples_parallel_8 = samples;
    Centroids centroids_parallel_8 = centroids;
    CentroidAssignment centroid_assignment_parallel_8 = centroid_assignment;

    auto start = std::chrono::high_resolution_clock::now();
    int iterations_1 = calculate_cluster(samples_serial, centroids_serial, centroid_assignment_serial, 1);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_simple = end - start;
    start = std::chrono::high_resolution_clock::now();
    calculate_cluster(samples_parallel_2, centroids_parallel_2, centroid_assignment_parallel_2, 2);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_parallel_2 = end - start;
    start = std::chrono::high_resolution_clock::now();
    calculate_cluster(samples_parallel_4, centroids_parallel_4, centroid_assignment_parallel_4, 4);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_parallel_4 = end - start;
    start = std::chrono::high_resolution_clock::now();
    calculate_cluster(samples_parallel_8, centroids_parallel_8, centroid_assignment_parallel_8, 8);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_parallel_8 = end - start;

	std::cout << "Elapsed time serial: " << elapsed_simple.count() << " seconds.\n";
	std::cout << "Elapsed time parallel (2 threads): " << elapsed_parallel_2.count() << " seconds.\n";
	std::cout << "Elapsed time parallel (4 threads): " << elapsed_parallel_4.count() << " seconds.\n";
	std::cout << "Elapsed time parallel (8 threads): " << elapsed_parallel_8.count() << " seconds.\n";

    std::ofstream out_file("task2_execution_results", std::ios::app);
	out_file << "Samples;Dimensions;Centroids;Max_Iterations;Iterations;1Thread;2Threads;4Threads;8Threads\n";
    out_file << n_samples << ";" << n_dimensions << ";" << n_centroids << ";" << max_iterations << ";"
             << iterations_1 << ";"
             << elapsed_simple.count() << ";"
             << elapsed_parallel_2.count() << ";"
             << elapsed_parallel_4.count() << ";"
		     << elapsed_parallel_8.count() << "\n";

    return 0;
}

