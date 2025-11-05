
#include <iostream>
#include <vector>
#include <random>
#include <array>
#include <limits>
#include <cmath>

constexpr std::size_t n_samples = 10000;
constexpr std::size_t n_dimensions = 3;
constexpr std::size_t n_centroids = 10;
const std::pair<float, float> space = { 0.0f, 100.0f };

// Feste, kompilierbare Definition eines "n-elementigen" Samples

using Point = std::array<float, n_dimensions>;
using Samples = std::array<Point, n_samples>;
using Centroids = std::array<Point, n_centroids>;
using CentroidAssignment = std::array<int, n_samples>;

void generate_initial_state(Samples& samples, Centroids& centroids, CentroidAssignment centroid_assignment) {
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

	centroid_assignment.fill(-1);
}

void assign_samples_to_centroid(const Samples& samples, const Centroids& centroids, CentroidAssignment& centroid_assignment) {
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

void update_centroids(const Samples& samples, Centroids& centroids, const CentroidAssignment& centroid_assignment) {
    // temporäre Summen initialisieren
    Centroids new_centroids;
    for (auto& p : new_centroids) p.fill(0.0f);
    std::array<int, n_centroids> counts = {};

    // Summiere alle Samples nach zugewiesenem Zentroid
    for (std::size_t i = 0; i < n_samples; ++i) {
        int idx = centroid_assignment[i];
        if (idx < 0 || idx >= static_cast<int>(n_centroids)) continue; // nicht zugewiesen
        ++counts[idx];
        for (std::size_t d = 0; d < n_dimensions; ++d)
            new_centroids[idx][d] += samples[i][d];
    }

    // Teile durch Anzahl der zugewiesenen Samples -> Mittelwert
    for (std::size_t j = 0; j < n_centroids; ++j) {
        if (counts[j] > 0) {
            for (std::size_t d = 0; d < n_dimensions; ++d)
                centroids[j][d] = new_centroids[j][d] / static_cast<float>(counts[j]);
        }
        // else: kein Sample dem Zentroid zugewiesen -> Zentroid bleibt unverändert
    }
}


int main()
{
    Samples samples;
	Centroids centroids;
	CentroidAssignment centroid_assignment;


}

