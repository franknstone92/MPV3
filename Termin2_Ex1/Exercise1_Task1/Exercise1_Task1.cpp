#include <iostream>
#include <vector>
#include <random>
#include <chrono>

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_int_distribution<> dist(1, 10);

void test_run(int m, int n);

void calc(int* matrix, int* vector, int* result, int m, int const n, int num_threads) {
    
    std::cout << "Calculating " << m << " x " << n << " matrix with " << num_threads << "Threads..." << std::endl;
    #pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < m; i++) {
        double line_res = 0;
        for (int j = 0; j < n; j++) {

            double field_res = matrix[i * n + j] * vector[j];
            line_res += field_res;
        }
        result[i] = line_res;
    }
}

int main(int argc, char** argv)
{
    int m;
    int n;
    int num_threads = 16;

    if (argc == 2) {
        num_threads = std::atoi(argv[1]);
    }

    test_run(8000000, 8);
    test_run(8000, 8000);
    test_run(8, 8000000);

}

void test_run(int m, int n) {
    int* matrix = new int[m * n];
    int* vector = new int[n];
    int* result1 = new int[m];
    int* result2 = new int[m];
    int* result4 = new int[m];
    int* result8 = new int[m];
    int* result16 = new int[m];
    int* result32 = new int[m];

    for (int i = 0; i < m * n; i++) {
        matrix[i] = dist(gen);
    }
    for (int i = 0; i < n; i++) {
        vector[i] = dist(gen);
    }


    auto start = std::chrono::high_resolution_clock::now();
    calc(matrix, vector, result1, m, n, 1);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed1 = end - start;

    start = std::chrono::high_resolution_clock::now();
    calc(matrix, vector, result2, m, n, 2);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed2 = end - start;

    start = std::chrono::high_resolution_clock::now();
    calc(matrix, vector, result4, m, n, 4);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed4 = end - start;

    start = std::chrono::high_resolution_clock::now();
    calc(matrix, vector, result8, m, n, 8);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed8 = end - start;

    start = std::chrono::high_resolution_clock::now();
    calc(matrix, vector, result16, m, n, 16);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed16 = end - start;

    start = std::chrono::high_resolution_clock::now();
    calc(matrix, vector, result32, m, n, 32);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed32 = end - start;

    std::cout << "Execution times for " << m << " x " << n << " calculation: " << std::endl;
    std::cout << " 1 core: " << elapsed1.count() << " ms" << std::endl;
    std::cout << " 2 core: " << elapsed2.count() << " ms" << std::endl;
    std::cout << " 4 core: " << elapsed4.count() << " ms" << std::endl;
    std::cout << " 8 core: " << elapsed8.count() << " ms" << std::endl;
    std::cout << "16 core: " << elapsed16.count() << " ms" << std::endl;
    std::cout << "32 core: " << elapsed32.count() << " ms" << std::endl;

    delete[] matrix;
    delete[] vector;
    delete[] result1;
    delete[] result2;
    delete[] result4;
    delete[] result8;
    delete[] result16;
    delete[] result32;
}


