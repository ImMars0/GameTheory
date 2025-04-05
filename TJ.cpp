#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <limits>
#include <algorithm>
using namespace std;

// Funcție pentru citirea matricei de plăți Q
void citesteMatriceQ(vector<vector<double>>& Q) {
    int randuri, coloane;
    cout << "Introdu numarul de strategii pentru jucatorul A (randuri): ";
    cin >> randuri;
    cout << "Introdu numarul de strategii pentru jucatorul B (coloane): ";
    cin >> coloane;

    Q.resize(randuri, vector<double>(coloane));
    cout << "Introdu elementele matricei Q :\n";
    for (int i = 0; i < randuri; ++i) {
        for (int j = 0; j < coloane; ++j) {
            cout << "Q[" << i+1 << "][" << j+1 << "] = ";
            cin >> Q[i][j];
        }
    }
}

// Funcție pentru afișarea matricei de plăți Q cu α și β
void afiseazaMatriceQ(const vector<vector<double>>& Q) {
    if (Q.empty() || Q[0].empty()) {
        cout << "Matricea este goala!\n";
        return;
    }

    int randuri = Q.size();
    int coloane = Q[0].size();

    // Afișăm antetul coloanelor
    cout << "\nMatricea Q cu alpha si betha:\n";
    cout << setw(8) << "|";
    for (int j = 0; j < coloane; ++j) {
        cout << setw(8) << "b_" << j+1 << " |";
    }
    cout << setw(8) << "alpha_i" << " |";
    cout << "\n" << string(8 + 12*(coloane+1), '-') << "\n";

    // Calculăm α_i (minimul pe fiecare rând)
    vector<double> alpha(randuri);
    for (int i = 0; i < randuri; ++i) {
        alpha[i] = *min_element(Q[i].begin(), Q[i].end());
    }

    // Afișăm fiecare rând cu valoarea α
    for (int i = 0; i < randuri; ++i) {
        cout << setw(6) << "a_" << i+1 << " |";
        for (int j = 0; j < coloane; ++j) {
            cout << setw(8) << Q[i][j] << " |";
        }
        cout << setw(8) << alpha[i] << " |";
        cout << "\n" << string(8 + 12*(coloane+1), '----') << "\n";
    }

    // Calculăm β_j (maximul pe fiecare coloană)
    vector<double> beta(coloane, -numeric_limits<double>::infinity());
    cout << setw(8) << " betha_j |";
    for (int j = 0; j < coloane; ++j) {
        for (int i = 0; i < randuri; ++i) {
            if (Q[i][j] > beta[j]) {
                beta[j] = Q[i][j];
            }
        }
        cout << setw(8) << beta[j] << " |";
    }
    cout << setw(8) << " |"; // Spațiu gol pentru coloana α
    cout << "\n" << string(8 + 12*(coloane+1), '----') << "\n";

    // Calculăm v1 și v2
    double v1 = *max_element(alpha.begin(), alpha.end());
    double v2 = *min_element(beta.begin(), beta.end());

    cout << "\nValoarea inferioara a jocului (v1): " << v1 << endl;
    cout << "Valoarea superioara a jocului (v2): " << v2 << endl;

    if (v1 == v2) {
        cout << "Jocul are punct de echilibru în strategii pure." << endl;
    } else {
        cout << "Jocul nu are punct de echilibru in strategii pure." << endl;
        cout << "Se cauta strategii mixte optime..." << endl;
    }
}

// Transformă problema jocului în problemă de programare liniară pentru jucătorul B (maximizare)
void transformaInPL(const vector<vector<double>>& Q, vector<vector<double>>& A, vector<double>& b, vector<double>& c) {
    int randuri = Q.size();
    int coloane = Q[0].size();

    // Jucătorul B maximizează v, cu restricțiile: Q^T * y ≤ 1, y ≥ 0
    // Funcția obiectiv: maximizează y_1 + y_2 + ... + y_n (echivalent cu maximizează v)
    A.resize(randuri, vector<double>(coloane, 0));
    b.resize(randuri, 1);
    c.resize(coloane, 1);

    for (int i = 0; i < randuri; ++i) {
        for (int j = 0; j < coloane; ++j) {
            A[i][j] = Q[i][j];
        }
    }
}


// Funcție pentru afișarea problemei de programare liniară pentru ambii jucatori
void afiseazaProblemaPL(const vector<vector<double>>& Q, bool pentruJucatorulB) {
    if (Q.empty() || Q[0].empty()) {
        cout << "Matricea este goala!\n";
        return;
    }

    int randuri = Q.size();
    int coloane = Q[0].size();

    if (pentruJucatorulB) {
        // Problema pentru jucătorul B (maximizare)
        cout << "\nProblema de programare liniara pentru JUCATORUL B (maximizare):\n";
       // cout << "Maximizeaza: v\n";
        cout << "Sub restrictiile:\n";
        
        for (int i = 0; i < randuri; ++i) {
            for (int j = 0; j < coloane; ++j) {
                cout << Q[i][j] << "y" << j+1;
                if (j < coloane-1) cout << " + ";
            }
            cout << " <= 1\n";
        }
        
        for (int j = 0; j < coloane; ++j) {
            cout << "y" << j+1 << " >= 0";
            if (j < coloane-1) cout << ", ";
        }
        cout << "\n";
    } else {
        // Problema pentru jucătorul A (minimizare)
        cout << "\nProblema de programare liniara pentru JUCATORUL A (minimizare):\n";
       // cout << "Minimizeaza: v\n";
        cout << "Sub restrictiile:\n";
        
        for (int j = 0; j < coloane; ++j) {
            for (int i = 0; i < randuri; ++i) {
                cout << Q[i][j] << "x" << i+1;
                if (i < randuri-1) cout << " + ";
            }
            cout << " >= 1\n";
        }
        
        for (int i = 0; i < randuri; ++i) {
            cout << "x" << i+1 << " >= 0";
            if (i < randuri-1) cout << ", ";
        }
        cout << "\n";
    }
}


// Aduce problema la forma standard pentru simplex
void formaStandard(vector<vector<double>>& A, vector<double>& b, vector<double>& c, 
                  vector<vector<double>>& A_std, vector<double>& b_std, vector<double>& c_std, 
                  vector<int>& baza) {
    int n = A.size();
    int m = A[0].size();

    // Adăugăm variabile de abatere (s_i) pentru fiecare restricție
    A_std.resize(n, vector<double>(m + n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            A_std[i][j] = A[i][j];
        }
        A_std[i][m + i] = 1;  // Coeficientul variabilei de abatere
    }

    b_std = b;
    c_std.resize(m + n, 0);
    for (int j = 0; j < m; ++j) {
        c_std[j] = c[j];
    }

    // Baza inițială conține variabilele de abatere
    baza.resize(n);
    for (int i = 0; i < n; ++i) {
        baza[i] = m + i;
    }
}

//Afisarea tabelului simplex inittial
// Functia pentru crearea tabelului simplex initial
vector<vector<double>> creare_tabel_simplex(
    const vector<vector<double>> & A_standard,
    const vector<double> & b_standard,
    const vector<double> & c_standard,
    const vector<int> & baza
) {
    int n = A_standard.size(); // numarul de restrictii
    int m = A_standard[0].size(); // numarul total de variabile
   
        

    // Tabelul simplex are dimensiunea (n+1) x (m+1)
    vector<vector<double>> tabel(n + 1, vector<double>(m + 1, 0));

    // Completam tabelul cu coeficientii restrictiilor si termenii liberi
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            tabel[i][j] = A_standard[i][j];
        }
        tabel[i][m] = b_standard[i]; // termenii liberi pe ultima coloana
    }

    // Calculam valorile Z_j si c_j - Z_j pentru fiecare coloana
    for (int j = 0; j < m; j++) {
        double z_j = 0;
        for (int i = 0; i < n; i++) {
            int var_baza = baza[i];
            z_j += c_standard[var_baza] * tabel[i][j];
        }
        tabel[n][j] = c_standard[j] - z_j; // c_j - z_j
    }

    // Calculam valoarea functiei obiectiv Z
    double Z = 0;
    for (int i = 0; i < n; i++) {
        int var_baza = baza[i];
        Z += c_standard[var_baza] * tabel[i][m];
    }
    tabel[n][m] = Z; // valoarea functiei obiectiv

    return tabel;
}

void afisare_tabel_simplex(
    const vector<vector<double>> & tabel,
    const vector<int> & baza,
    int iteratie,
    const vector<double> & c,
    pair<int, int> pivot = {-1, -1}  
) {
    int n = tabel.size() - 1; // numarul de restrictii
    int m = tabel[0].size() - 1; // numarul de variabile + 1

    cout << "\n\nTabelul simplex - Iteratia " << iteratie << ":\n";

    // Afisam antetul tabelului
    cout << setw(8) << "C_B" << " | " << setw(8) << "Baza" << " | " << setw(8) << "Solutia" << " | ";
    for (int j = 0; j < m; j++) {
        cout << setw(8) << "a" << j + 1 << " | ";
    }
    cout << endl;

    cout << string(10 + (m + 3) * 12, '-') << endl;

    // Afisam liniile restrictiilor
    for (int i = 0; i < n; i++) {
        cout << setw(8) << c[baza[i]] << " | " << setw(8) << "a" << baza[i] + 1 << " |  ";
        cout << setw(8) << fixed << setprecision(2) << tabel[i][m] << " | ";
        for (int j = 0; j < m; j++) {
            // Highlight the pivot element
            if (i == pivot.first && j == pivot.second) {
                cout << "[" << setw(6) << fixed << setprecision(2) << tabel[i][j] << "] |  ";
            } else {
                cout << setw(8) << fixed << setprecision(2) << tabel[i][j] << " |  ";
            }
        }
        cout << endl;
    }

    cout << string(10 + (m + 3) * 12, '-') << endl;

    // Afisam linia functiei obiectiv (c_j - z_j)
    cout << setw(8) << " " << " | " << setw(8) << " " << " | " << setw(8) << " " << " | ";
    for (int j = 0; j < m; j++) {
        cout << setw(8) << fixed << setprecision(2) << tabel[n][j] << "  |  ";
    }

    // Afisam valoarea functiei obiectiv Z
    double Z = tabel[n][m];
    cout <<"z_"<<iteratie<<"="<<setw(8) << fixed << setprecision(2) << Z << endl;

    // Afisam pivotul mai vizibil
    if (pivot.first != -1 && pivot.second != -1) {
        cout << "\n*** PIVOT: Elementul de pe linia " << pivot.first + 1 
             << ", coloana " << pivot.second + 1 
             << " (valoare: " << tabel[pivot.first][pivot.second] << ") ***" << endl;
    }
}


// Algoritmul simplex pentru maximizare cu afisare la fiecare iteratie
void simplex_cu_afisare(vector<vector<double>>& A, vector<double>& b, vector<double>& c, 
    vector<int>& baza, double& valoare, vector<double>& solutie,vector<double>& delta_j_absolut) {
int n = A.size();//numarul de restrictii
int m = A[0].size();//numarul de variabile
bool optim = false;
int iteratie = 0;
int m_b = m - n; //nr de coloane dat strategiile jucatorului B




// Afisam tabelul initial
vector<vector<double>> tabel = creare_tabel_simplex(A, b, c, baza);
//afisare_tabel_simplex(tabel, baza, iteratie, c);
//iteratie++;


while (!optim) {
    // Find entering variable
    int intrare = -1;
    double max_cj = -1e-9;
    vector<double> z_j(m, 0);
    
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
            z_j[j] += c[baza[i]] * A[i][j];
        }
        double cj_zj = c[j] - z_j[j];
        if (cj_zj > max_cj) {
            max_cj = cj_zj;
            intrare = j;
        }
    }
    
    if (max_cj <= 1e-9) {
        optim = true;
        break;
    }
    
    // Find leaving variable
    int iesire = -1;
    double min_raport = numeric_limits<double>::max();
    
    for (int i = 0; i < n; ++i) {
        if (A[i][intrare] > 0) {
            double raport = b[i] / A[i][intrare];
            if (raport < min_raport) {
                min_raport = raport;
                iesire = i;
            }
        }
    }
    
    if (iesire == -1) {
        cout << "Problema este nemărginită!\n";
        return;
    }
    
    // Display the current tableau with the pivot that will be used
    tabel = creare_tabel_simplex(A, b, c, baza);
    afisare_tabel_simplex(tabel, baza, iteratie, c, {iesire, intrare});
    iteratie++;
    
    // Store pivot element before modifying the matrix
    double pivot = A[iesire][intrare];
    
    // Update the basis
    baza[iesire] = intrare;
    
    // Perform pivoting
    for (int j = 0; j < m; ++j) {
        A[iesire][j] /= pivot;
    }
    b[iesire] /= pivot;
    
    for (int i = 0; i < n; ++i) {
        if (i != iesire) {
            double factor = A[i][intrare];
            for (int j = 0; j < m; ++j) {
                A[i][j] -= factor * A[iesire][j];
            }
            b[i] -= factor * b[iesire];
        }
    }
}

// Display the final tableau without pivot information
tabel = creare_tabel_simplex(A, b, c, baza);
afisare_tabel_simplex(tabel, baza, iteratie, c);

// Extragem soluția finală
valoare = 0;
solutie.resize(m, 0);
for (int i = 0; i < n; ++i) {
if (baza[i] < m) {
solutie[baza[i]] = b[i];
valoare += c[baza[i]] * b[i];
}
}

// Afisam tabelul final
tabel = creare_tabel_simplex(A, b, c, baza);
//afisare_tabel_simplex(tabel, baza, iteratie, c);

 // Calculam delta_j_absolut din tabelul final
 tabel = creare_tabel_simplex(A, b, c, baza);
 delta_j_absolut.resize(m);
 for (int j = 0; j < m; ++j) {
     delta_j_absolut[j] = abs(tabel[n][j]); // Linia n contine c_j - z_j
 }

 //calculam valoarea jocului
 double v = 1.0 / valoare;
 
cout << "\nSOLUTIE OPTIMA GASITA!\n";
cout << "\n=== REZULTATE ===\n\n";
cout << "Valoarea jocului: " << v << endl<<endl;
// In the simplex_cu_afisare function, replace the problematic part with:
cout << "Strategia optima pentru jucatorul A: \n";
int num_strategies_A = A.size(); // Number of rows in Q matrix (player A strategies)

// The slack variables in the simplex table correspond to player A's strategies
for(int i=0; i<num_strategies_A; i++) {
    // The slack variables start at index m_b in the delta_j_absolut vector
    if(m_b + i < m) {
        cout << "x" << (i+1) << " = " << delta_j_absolut[m_b + i] * v << " (" 
             << round(delta_j_absolut[m_b + i] * v * 100) << "%)\n";
    }
}



}

int main() {
    vector<vector<double>> Q;
    citesteMatriceQ(Q);
    afiseazaMatriceQ(Q);


    // Transformăm problema pentru jucătorul B (maximizare)
    vector<vector<double>> A;
    vector<double> b, c;
    transformaInPL(Q, A, b, c);

    //Afi problemele PL asociate celor  două jucători
    afiseazaProblemaPL(Q,true);//B
    afiseazaProblemaPL(Q,false);//A

cout<<"\n\n NE VOM OCUPA DE PROBLEMA DE MAXIMIZARE CORESPUNZATOARE JUCATORULUI B"<<endl;

    // Aducem la forma standard
    vector<vector<double>> A_std;
    vector<double> b_std, c_std;
    vector<int> baza;
    formaStandard(A, b, c, A_std, b_std, c_std, baza);

    // Rezolvăm cu simplex
    double valoare;
    vector<double> solutie;
    vector<double> delta_final;
    simplex_cu_afisare(A_std, b_std, c_std, baza, valoare, solutie, delta_final);
//afisam delta_j_absolut

for (int j = Q[0].size(); j < A[0].size(); ++j) {
    cout << "delta_" << j + 1 << " = " << delta_final[j] << endl;
}
    // Calculăm strategiile mixte și valoarea jocului
    double v = 1.0 / valoare;
    vector<double> y(solutie.begin(), solutie.begin() + Q[0].size());
    for (double& yi : y) yi *= v;
    

    // Afișăm rezultatele
  
    cout << "\nStrategia optima pentru jucatorul B:\n";
    for (int j = 0; j < y.size(); ++j) {
        cout << "y" << j+1 << " = " << y[j] << " (" << round(y[j]*100) << "%)\n";
    }
   

    return 0;
}