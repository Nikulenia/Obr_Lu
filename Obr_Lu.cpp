#include <iostream>
#include <vector>
using namespace std;


// Функция для нахождения LU-разложения матрицы
void LU_decomposition(vector<vector<double>>& A, vector< vector<double>>& L, vector<vector<double>>& U) {
    int n = A.size();
    for (int i = 0; i < n; i++) {
        for (int k = i; k < n; k++) {
            double sum = 0;
            for (int j = 0; j < i; j++) {
                sum += (L[i][j] * U[j][k]);
            }
            U[i][k] = A[i][k] - sum;
        }
        for (int k = i; k < n; k++) {
            if (i == k) {
                L[i][i] = 1;
            }
            else {
                double sum = 0;
                for (int j = 0; j < i; j++) {
                    sum += (L[k][j] * U[j][i]);
                }
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }
}

// Функция для решения системы уравнений
vector<double> solve_LU(vector<vector<double>>& L, vector<vector<double>>& U, vector<int>& P, vector<double>& b) {
    int n = L.size();
    vector<double> y(n, 0);
    vector<double> x(n, 0);

    // Решаем Ly = P*b
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b[P[i]] - sum) / L[i][i];
    }

    // Решаем Ux = y
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }

    return x;
}

// Функция для нахождения обратной матрицы
vector<vector<double>> inverse_matrix(vector<vector<double>>& A) {
    int n = A.size();
    vector<vector<double>> L(n, vector<double>(n, 0));
    vector<vector<double>> U(n, vector<double>(n, 0));
    vector<int> P(n);
    for (int i = 0; i < n; i++) {
        P[i] = i;
    }

    LU_decomposition(A, L, U);

    vector<vector<double>> A_inv(n, vector<double>(n, 0));

    for (int i = 0; i < n; i++) {
        vector<double> b(n, 0);
        b[i] = 1;

        vector<double> x = solve_LU(L, U, P, b);

        for (int j = 0; j < n; j++) {
            A_inv[j][i] = x[j];
        }
    }

    return A_inv;
}

int main() {
    setlocale(LC_ALL, "ru");

    int n;
    cout << "Введите размерность матрицы: ";
    cin >> n;

    cout << "Введите элементы матрицы A: " << endl;
    vector<vector<double>> A(n, vector<double>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cin >> A[i][j];
        }
    }

    vector<vector<double>> A_inv = inverse_matrix(A);

    cout << "Обратная матрица:\t" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << A_inv[i][j] << " ";
        }
        cout << endl;
    }

    return 0;
}
