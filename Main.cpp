#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

long long getGcd(long long a, long long b) {
    a = abs(a); b = abs(b);
    while (b) {
        a %= b;
        swap(a, b);
    }
    return a;
}

struct Frac {
    long long num, den;
    
    Frac(long long n = 0, long long d = 1) {
        long long g = getGcd(n, d);
        num = n / g;
        den = d / g;
        if (den < 0) { num = -num; den = -den; }
    }
    
    Frac operator+(const Frac& o) const { return Frac(num * o.den + o.num * den, den * o.den); }
    Frac operator-(const Frac& o) const { return Frac(num * o.den - o.num * den, den * o.den); }
    Frac operator*(const Frac& o) const { return Frac(num * o.num, den * o.den); }
    Frac operator/(const Frac& o) const { return Frac(num * o.den, den * o.num); }
    bool operator==(const Frac& o) const { return num == o.num && den == o.den; }
    bool operator!=(const Frac& o) const { return !(*this == o); }
};

ostream& operator<<(ostream& os, const Frac& f) {
    if (f.den == 1 || f.num == 0) os << f.num;
    else os << f.num << "/" << f.den;
    return os;
}

typedef vector<vector<Frac>> Mat;

Mat getMinor(Mat m, int r, int c) {
    Mat res;
    for (int i = 0; i < m.size(); i++) {
        if (i == r) continue;
        vector<Frac> row;
        for (int j = 0; j < m.size(); j++) {
            if (j == c) continue;
            row.push_back(m[i][j]);
        }
        res.push_back(row);
    }
    return res;
}

Frac determinant(Mat m) {
    if (m.size() == 1) return m[0][0];
    if (m.size() == 2) return m[0][0] * m[1][1] - m[0][1] * m[1][0];
    Frac det(0, 1);
    for (int c = 0; c < m.size(); c++) {
        Frac term = m[0][c] * determinant(getMinor(m, 0, c));
        if (c % 2 == 1) det = det - term;
        else det = det + term;
    }
    return det;
}

Mat inverseDet(Mat m, Frac det) {
    int n = m.size();
    Mat adj(n, vector<Frac>(n));
    if (n == 1) {
        adj[0][0] = Frac(1, 1) / det;
        return adj;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            Frac sign = ((i + j) % 2 == 0) ? Frac(1, 1) : Frac(-1, 1);
            adj[j][i] = (determinant(getMinor(m, i, j)) * sign) / det;
        }
    }
    return adj;
}

Mat inverseGJ(Mat m) {
    int n = m.size();
    Mat aug(n, vector<Frac>(2 * n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) aug[i][j] = m[i][j];
        aug[i][n + i] = Frac(1, 1);
    }
    
    for (int i = 0; i < n; i++) {
        if (aug[i][i].num == 0) {
            bool swapped = false;
            for (int k = i + 1; k < n; k++) {
                if (aug[k][i].num != 0) {
                    swap(aug[i], aug[k]);
                    swapped = true;
                    break;
                }
            }
            if (!swapped) return Mat();
        }
        Frac p = aug[i][i];
        for (int j = 0; j < 2 * n; j++) aug[i][j] = aug[i][j] / p;
        for (int k = 0; k < n; k++) {
            if (k != i) {
                Frac f = aug[k][i];
                for (int j = 0; j < 2 * n; j++) {
                    aug[k][j] = aug[k][j] - (f * aug[i][j]);
                }
            }
        }
    }
    
    Mat inv(n, vector<Frac>(n));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i][j] = aug[i][n + j];
    return inv;
}

bool compareMat(Mat m1, Mat m2) {
    if (m1.empty() || m2.empty()) return m1.empty() == m2.empty();
    for (int i = 0; i < m1.size(); i++)
        for (int j = 0; j < m1.size(); j++)
            if (m1[i][j] != m2[i][j]) return false;
    return true;
}

void printMat(Mat m) {
    for (int i = 0; i < m.size(); i++) {
        for (int j = 0; j < m.size(); j++) {
            cout << m[i][j] << " ";
        }
        cout << "\n";
    }
}

int main() {
    int n;
    cout << "정수 n 입력: ";
    if (!(cin >> n)) return 0;
    cout << n << "x" << n << " 행렬의 원소를 행 단위로 공백으로 구분하여 입력하세요:\n";
    
    Mat m(n, vector<Frac>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            long long val;
            cin >> val;
            m[i][j] = Frac(val, 1);
        }
    }
    
    Frac det = determinant(m);
    if (det.num == 0) {
        cout << "\n오류: 역행렬이 존재하지 않습니다 (행렬식이 0입니다).\n";
        return 0;
    }
    
    Mat inv1 = inverseDet(m, det);
    cout << "\n[방법 1] 행렬식을 이용한 역행렬 (분수 표현):\n";
    printMat(inv1);
    
    Mat inv2 = inverseGJ(m);
    cout << "\n[방법 2] 가우스-조던 소거법을 이용한 역행렬 (분수 표현):\n";
    printMat(inv2);
    
    cout << "\n[결과 비교]\n";
    if (compareMat(inv1, inv2)) cout << "두 방법으로 계산한 역행렬이 동일합니다.\n";
    else cout << "두 방법의 결과가 다릅니다.\n";
    
    return 0;
}