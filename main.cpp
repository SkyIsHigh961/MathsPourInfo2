#include <iostream>
#include <string>
#include <vector>   
#include <cmath> //for additional math functions
#include <stdexcept> //for std::invalid_argument 
#include <fstream> //for file input/output
#include <sstream> //for std::stringstream

//******************************************************************************************************************************

using namespace std;

//******************************************************************************************************************************

vector<long double> addition(vector<long double> v1, vector<long double> v2){

    //variables
    vector<long double> result(v1.size());

    //vérfication de la taille des vecteurs
    if(v1.size() != v2.size()){
        throw invalid_argument("Les vecteurs ne sont pas de même taille");
    }

    //addition des vecteurs
    for (size_t i = 0; i < v1.size(); i++)
    {
        result[i]= v1[i] + v2[i];
    }
    
    return result;
}

//******************************************************************************************************************************

vector<long double> multiplication_scalaire(vector<long double> v, long double scalar) {
    
    //variables
    vector<long double> result(v.size());

    // Multiplication de chaque élément du vecteur par le scalaire
    for (size_t i = 0; i < v.size(); i++) {
        result[i] = v[i] * scalar;
    }

    return result;
}

//******************************************************************************************************************************

 vector<long double> produit_scalaire(vector<long double> v1, vector<long double> v2){
    
    //variables
    vector<long double> result(v1.size());

    //vérfication de la taille des vecteurs
    if(v1.size() != v2.size()){
        throw invalid_argument("Les vecteurs ne sont pas de même taille");
    }

    //scalar product of the vectors
    for (size_t i = 0; i < v1.size(); i++)
    {
        result[i]= v1[i] * v2[i];
    }
    
    return result;
 }

//******************************************************************************************************************************
// Fonction pour transposer une matrice

vector<vector<long double>> transpose_matrice(const vector<vector<long double>>& matrix) {
    if (matrix.empty()) return {};

    // Crée une nouvelle matrice avec les dimensions inversées
    vector<vector<long double>> transposed(matrix[0].size(), vector<long double>(matrix.size()));

    // Remplissage de la matrice transposée
    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix[0].size(); j++) {
            transposed[j][i] = matrix[i][j];
        }
    }

    return transposed;
}

//******************************************************************************************************************************

vector<vector<long double>> produit_matrices(const vector<vector<long double>>& mat1, const vector<vector<long double>>& mat2) {
    
    // Vérification de la compatibilité des matrices pour la multiplication
    if (mat1.empty() || mat2.empty() || mat1[0].size() != mat2.size()) {
        throw invalid_argument("Les matrices ne peuvent pas être multipliées en raison de tailles incompatibles");
    }

    // Initialisation de la matrice résultante avec des zéros
    vector<vector<long double>> result(mat1.size(), vector<long double>(mat2[0].size(), 0));

    // Boucle sur chaque ligne de mat1
    for (size_t i = 0; i < mat1.size(); ++i) {
        // Boucle sur chaque colonne de mat2
        for (size_t j = 0; j < mat2[0].size(); ++j) {
            // Calcul de l'élément [i][j] du résultat
            for (size_t k = 0; k < mat2.size(); ++k) {
                result[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }

    return result;
}

//******************************************************************************************************************************

vector<long double> appliquer_matrice_vecteur(const vector<vector<long double>>& mat, const vector<long double>& vec) {
    
    // Vérifier si la matrice est vide ou si le nombre de colonnes de la matrice ne correspond pas à la taille du vecteur
    if (mat.empty() || mat[0].size() != vec.size()) {
        throw invalid_argument("Taille incompatibles entre la matrice et le vecteur");
    }

    // Initialisation du vecteur résultant
    vector<long double> result(mat.size(), 0.0);

    // Boucle sur chaque ligne de la matrice
    for (size_t i = 0; i < mat.size(); i++) {
        // Boucle sur chaque colonne de la ligne actuelle
        for (size_t j = 0; j < mat[i].size(); j++) {
            // Multiplication de l'élément de la matrice par l'élément correspondant du vecteur et addition au résultat
            result[i] += mat[i][j] * vec[j];
        }
    }

    return result;
}

//******************************************************************************************************************************

bool isValide(vector<vector<long double>> matrice){
    for (vector<long double> e: matrice){
        if(matrice[0].size() != e.size()){
            std::cout<< "matrice invalide"<< std::endl;
            return 0;
        }
    }
    std::cout<< "matrice valide"<< std::endl;
    return 1;
}

//******************************************************************************************************************************
// Cette fonction réalise la décomposition de Cholesky d'une matrice symétrique définie positive
// Précondition : La matrice en entrée doit etre définie positive sinon la fonction renvoie une erreur.

vector<vector<long double>> cholesky_decomposition(const vector<vector<long double>>& A) {
    size_t n = A.size();
    vector<vector<long double>> B(n, vector<long double>(n, 0.0));

    for (size_t j = 0; j < n; j++) {
        long double sum = 0.0;
        for (size_t k = 0; k < j; k++) {
            sum += B[j][k] * B[j][k];
        }
        long double diag = A[j][j] - sum;
        if (diag <= 0) {
            throw std::runtime_error("La matrice n'est pas définie positive.");
        }
        B[j][j] = sqrt(diag);

        for (size_t i = j + 1; i < n; i++) {
            sum = 0.0;
            for (size_t k = 0; k < j; k++) {
                sum += B[i][k] * B[j][k];
            }
            B[i][j] = (A[i][j] - sum) / B[j][j];
        }
    }

    return B;
}

//******************************************************************************************************************************

vector<long double> substitution_avant(const vector<vector<long double>>& L, const vector<long double>& b) {
    size_t n = L.size();
    vector<long double> y(n, 0.0);

    for (size_t i = 0; i < n; ++i) {
        y[i] = b[i];
        for (size_t j = 0; j < i; ++j) {
            y[i] -= L[i][j] * y[j];
        }
        y[i] /= L[i][i];
    }

    return y;
}

//******************************************************************************************************************************

vector<long double> substitution_arriere(const vector<vector<long double>>& U, const vector<long double>& y) {
    size_t n = U.size();
    vector<long double> x(n, 0.0);

    for (int i = n - 1; i >= 0; --i) {
        x[i] = y[i];
        for (size_t j = i + 1; j < n; ++j) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }

    return x;
}

//******************************************************************************************************************************
// Function to calculate the solution of the normal equation A^TAx = A^Tb

vector<long double> resoudre_equation_normale(const vector<vector<long double>>& A, const vector<long double>& b) {
    
    // Vérifiez que le nombre de lignes de A correspond à la taille de b
    if (A.empty() || A.size() != b.size()) {
        throw invalid_argument("La matrice et le vecteur ne sont pas compatibles pour résoudre l'équation normale.");
    }
    
    vector<vector<long double>> At = transpose_matrice(A); // At est la transposée de A
    vector<vector<long double>> AtA = produit_matrices(At, A); // AtA est le produit de At par A
    vector<long double> Atb = appliquer_matrice_vecteur(At, b); // Atb est le produit de At par b

    vector<vector<long double>> B = cholesky_decomposition(AtA); // B est la décomposition de Cholesky de AtA
    vector<vector<long double>> Bt = transpose_matrice(B); // Bt est la transposée de B

    // Résoudre By = Atb pour y par substitution avant
    vector<long double> y = substitution_avant(B, Atb);
    // Résoudre B^Tx = y pour x par substitution arrière
    vector<long double> x_chapeau = substitution_arriere(Bt, y);

    return x_chapeau;
}

//Test :
// Définissons une matrice A et un vecteur b
    // vector<vector<long double>> A = {
    //     {2, -1, 0, 0},
    //     {-1, 2, -1, 0},
    //     {0, -1, 2, -1},
    //     {0, 0, -1, 2}
    // };
    // vector<long double> b = {1, 0, 0, 0};

    // // Utilisons la fonction pour calculer la solution des moindres carrés
    // try {
    //     vector<long double> x_chapeau = resoudre_equation_normale(A, b);
        
    //     // Affichons le résultat
    //     cout << "La solution x_chapeau est:" << endl;
    //     for (long double val : x_chapeau) {
    //         cout << val << " ";
    //     }
    //     cout << endl;
    // } catch (const std::exception& e) {
    //     cerr << "Erreur: " << e.what() << endl;
    // }

    // return 0;

//******************************************************************************************************************************

void readHousingData(const string& filename, vector<vector<long double>>& A, vector<long double>& b,
                      bool includeArea = false, bool includeBedrooms = false,
                      bool includeBathrooms = false, bool includeCondoStatus = false,
                      int propertyTypeFilter = -1) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Could not open file: " + filename);
    }

    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        string token;
        vector<string> tokens;

        while (getline(ss, token, ',')) {
            tokens.push_back(token);
        }

        if (tokens.size() < 6) continue;  // Ensure there are enough columns

        // Property type check
        int propertyType = stoi(tokens[4]);  // Index 4 for property type (-1 = aucun, 0 = maison individuelles, 1 = copropriétés)
        if (propertyTypeFilter != -1 && propertyType != propertyTypeFilter) continue;

        vector<long double> rowA{1.0};  // Start with the intercept term

        if (includeArea) {
            rowA.push_back(stod(tokens[3]));  // Include surface area (index 3)
        }

        if (includeBedrooms) {
            rowA.push_back(stod(tokens[1]));  // Include number of bedrooms (index 1)
        }

        if (includeBathrooms) {
            rowA.push_back(stod(tokens[2]));  // Include number of bathrooms (index 2)
        }

        if (includeCondoStatus) {
            rowA.push_back(stod(tokens[4]));  // Include condo status (index 4, 0 or 1)
        }

        A.push_back(rowA);
        b.push_back(stod(tokens[5]));  // Price is always included (index 5)
    }

    file.close();
}

    //For testing exercice 2.3
    //variables
    // vector<vector<long double>> A;
    // vector<long double> b;
    // string filename = "housing.data.txt";

    // readHousingData(filename, A, b);

    // vector<long double> x_chapeau = resoudre_equation_normale(A, b);

    //     // Affichons le résultat
    //     cout << "La solution x_chapeau est:" << endl;
    //     for (long double val : x_chapeau) {
    //         cout << val << " ";
    //     }
    //     cout << endl;


    // return 0;

//******************************************************************************************************************************

long double calculeSigmaChapeau(const vector<vector<long double>>& A, const vector<long double>& b, const vector<long double>& theta) {
    vector<long double> residus(b.size());
    long double sommeDesCarresDesResidus = 0.0;

    // Calculer les valeurs prédites
    for (size_t i = 0; i < A.size(); ++i) {
        long double predite = 0.0;
        for (size_t j = 0; j < A[i].size(); ++j) {
            predite += A[i][j] * theta[j];
        }
        residus[i] = b[i] - predite;
        sommeDesCarresDesResidus += std::pow(residus[i], 2);
    }

    // Calculer sigma chapeau
    return std::sqrt(sommeDesCarresDesResidus / b.size());
}

    // vector<vector<long double>> A;
    // vector<long double> b;
    // string filename = "housing.data.txt"; // Ensure this path is correct

    // readHousingData(filename, A, b);
    // vector<long double> theta = resoudre_equation_normale(A, b);
    // long double sigmaHat = calculeSigmaChapeau(A, b, theta);

    // cout << "Sigma Chapeau: " << sigmaHat << endl;

    // return 0; 

//******************************************************************************************************************************
// For debug purposes
void writeMatrixToFile(const vector<vector<long double>>& matrix, const string& filename) {
    ofstream outFile(filename);  // Create an ofstream to write to the file
    if (!outFile.is_open()) {
        cerr << "Failed to open file: " << filename << endl;
        return;
    }

    // Loop through each row of the matrix
    for (const auto& row : matrix) {
        for (size_t i = 0; i < row.size(); i++) {
            outFile << row[i];  // Write each element to the file
            if (i != row.size() - 1) outFile << ", ";  // Separate elements with a comma
        }
        outFile << "\n";  // Start a new line after each row
    }

    outFile.close();  // Close the file after writing
    cout << "Matrix written to " << filename << endl;
}

//******************************************************************************************************************************\

// Function to generate polynomial terms for a vector of variables
void generatePolynomialTerms(const vector<long double>& vars, int degree, vector<long double>& row) {
    row.push_back(1.0);  // Add intercept term
    for (int i = 1; i <= degree; ++i) {
        for (long double var : vars) {
            row.push_back(pow(var, i));  // Add variable raised to the power 'i'
        }
    }
}

// Main function to read data and create the matrix A and vector b
void createPolynomialRegressionMatrix(vector<vector<long double>>& data, int degree, vector<vector<long double>>& A) {
    for (const auto& row : data) {
        vector<long double> A_row;
        vector<long double> vars(row.begin(), row.end() - 1);  // Extract all variables (assuming last column is the dependent variable)
        generatePolynomialTerms(vars, degree, A_row);
        A.push_back(A_row);
        //b.push_back(row.back());  // Last element is the dependent variable
    }
}


void readpolyData(const string& filename, vector<vector<long double>>& A, vector<long double>& b) {
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("Could not open file: " + filename);
    }

    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        string token;
        vector<long double> rowA;  // Create a vector for the current row of matrix A
        long double value;  // Variable to store the parsed long double value from tokens

        // Read each token separated by comma
        while (getline(ss, token, ',')) {
            // Convert the token to long double
            value = stod(token);
            // Add the value to the row of matrix A
            rowA.push_back(value);
        }

        // Check if the row contains at least three values (including the intercept term)
        if (rowA.size() < 3) {
            throw runtime_error("Insufficient columns in the data file.");
        }

        // Add the row to matrix A
        A.push_back(rowA);

        // Add the value from the third column to vector b
        b.push_back(rowA.back());
    }

    file.close();
}




int main(){
    
    vector<vector<long double>> A;
    vector<vector<long double>> newA;
    vector<long double> b;
    int polynomialDegree=3;
    //readHousingData("housing.data.txt", A, b, true, false, false, false, -1);  
    readpolyData("polynomial.data.txt",A,b);
    
    createPolynomialRegressionMatrix(A, polynomialDegree, newA);

    vector<long double> theta = resoudre_equation_normale(newA, b);
    long double sigmaHat = calculeSigmaChapeau(newA, b, theta);

    cout <<"pour d= "<<polynomialDegree<< " Sigma Chapeau:  " << sigmaHat << endl;
    
    
}