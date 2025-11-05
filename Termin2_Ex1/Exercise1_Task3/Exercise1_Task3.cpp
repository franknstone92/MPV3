#include <iostream>
#include <random>
#include <fstream>
#include <format>
#include <chrono>

struct Vector {
    double x;
    double y;
    double z;
};

struct Particle {
    double m;
    Vector s; // position
    Vector v; // velocity
};

const double G = 6.6743e-11;

const int THREADS_1 = 1;
const int THREADS_2 = 2;
const int THREADS_4 = 4;
const int THREADS_8 = 8;
const int THREADS_16 = 16;
const int THREADS_32 = 32;

int n_steps = 10000;
int n_particles = 200;
int delta_t = 1e4;  // hours

void generate_initial_state(Particle* particles) {
    std::default_random_engine generator{ 42 };

    std::lognormal_distribution<double> mass_distibution{ std::log(1e25), 0.8 };
    std::uniform_real_distribution<double> position_distribution{ -1e11, +1e11 };
    std::normal_distribution<double> velocity_distribution{ 0.0, 1e2 };
    
    for (int q = 0; q < n_particles; q++) {
        particles[q].m = mass_distibution(generator);
        
        particles[q].s.x = position_distribution(generator);
        particles[q].s.y = position_distribution(generator);
        particles[q].s.z = position_distribution(generator);
        
        particles[q].v.x = position_distribution(generator);
        particles[q].v.y = position_distribution(generator);
        particles[q].v.z = position_distribution(generator);
    }
}

// simple calculation
void calculate_forces_simple(int q,  Particle* particles, Vector* forces) {
    forces[q].x = forces[q].y = forces[q].z = 0;

    for (int k = 0; k < n_particles; k++) {
        if (k == q) continue;

        double dx = particles[q].s.x - particles[k].s.x;
        double dy = particles[q].s.y - particles[k].s.y;
        double dz = particles[q].s.z - particles[k].s.z;

        double r2 = dx * dx + dy * dy + dz * dz;
        double r1 = std::sqrt(r2);
        double r3 = r2 * r1;

        double mass_force = particles[k].m / r3;

        forces[q].x += mass_force * dx;
        forces[q].y += mass_force * dy;
        forces[q].z += mass_force * dz;
    }

    forces[q].x *= -G;
    forces[q].y *= -G;
    forces[q].z *= -G;
}

// reduced calculation using newton's third law
void calculate_forces_reduced(int q, Particle* particles, Vector* forces) {

    // compute pairwise forces only once
    for (int k = q + 1; k < n_particles; k++) {
        double dx = particles[q].s.x - particles[k].s.x;
        double dy = particles[q].s.y - particles[k].s.y;
        double dz = particles[q].s.z - particles[k].s.z;

        double r2 = dx * dx + dy * dy + dz * dz;
        double r1 = std::sqrt(r2);
        double r3 = r2 * r1;

        double mass_force = particles[k].m / r3;

        
        forces[q].x += mass_force * dx;
        forces[q].y += mass_force * dy;
        forces[q].z += mass_force * dz;

        forces[k].x -= mass_force * dx;
        forces[k].y -= mass_force * dy;
        forces[k].z -= mass_force * dz;
    }
}


void update_position(int q, Particle* particles, Vector* forces) {
    particles[q].s.x += particles[q].v.x * delta_t;
    particles[q].s.y += particles[q].v.y * delta_t;
    particles[q].s.z += particles[q].v.z * delta_t;

    particles[q].v.x += forces[q].x * delta_t;
    particles[q].v.y += forces[q].y * delta_t;
    particles[q].v.z += forces[q].z * delta_t;
}

void write_header(std::ostream& output) {
    output << "Step;Time;Particle;Position_X;Position_Y;Position_Z\n";

}

void write_state(std::ostream& output, int step, double time, Particle* particles) {
    for (int q = 0; q < n_particles; q++) {
        output << std::format("{};{};{};{};{}", step, time, q, particles[q].s.x, particles[q].s.y, particles[q].s.z);

    }
}

void simple_calculation(int num_threads) {
    Particle* particles = new Particle[n_particles];
    generate_initial_state(particles);

    std::fstream out_file{ "n_body_serial_simple.csv", std::ios::out };
    write_header(out_file);
    write_state(out_file, 0, 0.0, particles);

    Vector* forces = new Vector[n_particles];

    for (int step = 1; step < n_steps; step++) {
        double current_time = step * delta_t;

        // calculate forces simple
        #pragma omp parallel for num_threads(num_threads) //schedule(static, 10)
        for (int q = 0; q < n_particles; q++)
            calculate_forces_simple(q, particles, forces);

        // position updates
        for (int q = 0; q < n_particles; q++) {
            update_position(q, particles, forces);
        }

        // write results
        if (step % 100 == 0) {
            std::cout << std::format("{} / {}\n", step, n_steps);
            write_state(out_file, step, current_time, particles);
        }
    }

    delete[] particles;
}

void reduced_calculation(int num_threads) {
    Particle* particles = new Particle[n_particles];
    generate_initial_state(particles);

    std::fstream out_file{ "n_body_serial_reduced.csv", std::ios::out };
    write_header(out_file);
    write_state(out_file, 0, 0.0, particles);

    Vector* forces = new Vector[n_particles];
    std::cout << "Calculating using " << num_threads << " Threads." << std::endl;

    for (int step = 1; step < n_steps; step++) {
        double current_time = step * delta_t;

        for (int i = 0; i < n_particles; i++) {
            forces[i].x = forces[i].y = forces[i].z = 0;
        }

        // calculate forces reduced
        #pragma omp parallel for num_threads(num_threads) reduction(+:forces[:n_particles]) schedule(dynamic, 10)
        for (int q = 0; q < n_particles; q++)
            calculate_forces_reduced(q, particles, forces);

        // position updates
        for (int q = 0; q < n_particles; q++) {
            update_position(q, particles, forces);
        }

        // write results
        if (step % 100 == 0) {
            std::cout << std::format("{} / {}\n", step, n_steps);
            write_state(out_file, step, current_time, particles);
        }
    }

    delete[] particles;
}


int main()
{
    // Serial calculation
    auto start = std::chrono::high_resolution_clock::now();
    simple_calculation(THREADS_1);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_simple = end - start;
    
    start = std::chrono::high_resolution_clock::now();
    reduced_calculation(THREADS_1);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_reduced = end - start;

    // Parallel calculation
    start = std::chrono::high_resolution_clock::now();
    simple_calculation(THREADS_8);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_simple_8_threads = end - start;

    start = std::chrono::high_resolution_clock::now();
    reduced_calculation(THREADS_8);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_reduced_8_threads = end - start;

    std::cout << "Execution times for SERIAL " << n_particles << "-body calculation: " << std::endl;
    std::cout << " Simple calculation : " << elapsed_simple.count() << " seconds" << std::endl;
    std::cout << " Reduced calculation: " << elapsed_reduced.count() << " seconds" << std::endl;

    std::cout << "Execution times for PARALLEL " << n_particles << "-body calculation: " << std::endl;
    std::cout << " Simple calculation with  8 Threads: " << elapsed_simple_8_threads.count() << " seconds" << std::endl;
    std::cout << " Reduced calculation with 8 Threads: " << elapsed_reduced_8_threads.count() << " seconds" << std::endl;

}
