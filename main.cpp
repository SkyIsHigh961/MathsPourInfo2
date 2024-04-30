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

vector<double> addition(vector<double> v1, vector<double> v2){

    //variables
    vector<double> result(v1.size());

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

vector<double> multiplication_scalaire(vector<double> v, double scalar) {
    
    //variables
    vector<double> result(v.size());

    // Multiplication de chaque élément du vecteur par le scalaire
    for (size_t i = 0; i < v.size(); i++) {
        result[i] = v[i] * scalar;
    }

    return result;
}

//******************************************************************************************************************************

 vector<double> produit_scalaire(vector<double> v1, vector<double> v2){
    
    //variables
    vector<double> result(v1.size());

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

vector<vector<double>> transpose_matrice(const vector<vector<double>>& matrix) {
    if (matrix.empty()) return {};

    // Crée une nouvelle matrice avec les dimensions inversées
    vector<vector<double>> transposed(matrix[0].size(), vector<double>(matrix.size()));

    // Remplissage de la matrice transposée
    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix[0].size(); j++) {
            transposed[j][i] = matrix[i][j];
        }
    }

    return transposed;
}

//******************************************************************************************************************************

vector<vector<double>> produit_matrices(const vector<vector<double>>& mat1, const vector<vector<double>>& mat2) {
    
    // Vérification de la compatibilité des matrices pour la multiplication
    if (mat1.empty() || mat2.empty() || mat1[0].size() != mat2.size()) {
        throw invalid_argument("Les matrices ne peuvent pas être multipliées en raison de tailles incompatibles");
    }

    // Initialisation de la matrice résultante avec des zéros
    vector<vector<double>> result(mat1.size(), vector<double>(mat2[0].size(), 0));

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

vector<double> appliquer_matrice_vecteur(const vector<vector<double>>& mat, const vector<double>& vec) {
    
    // Vérifier si la matrice est vide ou si le nombre de colonnes de la matrice ne correspond pas à la taille du vecteur
    if (mat.empty() || mat[0].size() != vec.size()) {
        throw invalid_argument("Taille incompatibles entre la matrice et le vecteur");
    }

    // Initialisation du vecteur résultant
    vector<double> result(mat.size(), 0.0);

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

bool isValide(vector<vector<double>> matrice){
    for (vector<double> e: matrice){
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

vector<vector<double>> cholesky_decomposition(const vector<vector<double>>& A) {
    size_t n = A.size();
    vector<vector<double>> B(n, vector<double>(n, 0.0));

    for (size_t j = 0; j < n; j++) {
        double sum = 0.0;
        for (size_t k = 0; k < j; k++) {
            sum += B[j][k] * B[j][k];
        }
        double diag = A[j][j] - sum;
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

vector<double> substitution_avant(const vector<vector<double>>& L, const vector<double>& b) {
    size_t n = L.size();
    vector<double> y(n, 0.0);

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

vector<double> substitution_arriere(const vector<vector<double>>& U, const vector<double>& y) {
    size_t n = U.size();
    vector<double> x(n, 0.0);

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

vector<double> resoudre_equation_normale(const vector<vector<double>>& A, const vector<double>& b) {
    
    // Vérifiez que le nombre de lignes de A correspond à la taille de b
    if (A.empty() || A.size() != b.size()) {
        throw invalid_argument("La matrice et le vecteur ne sont pas compatibles pour résoudre l'équation normale.");
    }
    
    vector<vector<double>> At = transpose_matrice(A); // At est la transposée de A
    vector<vector<double>> AtA = produit_matrices(At, A); // AtA est le produit de At par A
    vector<double> Atb = appliquer_matrice_vecteur(At, b); // Atb est le produit de At par b

    vector<vector<double>> B = cholesky_decomposition(AtA); // B est la décomposition de Cholesky de AtA
    vector<vector<double>> Bt = transpose_matrice(B); // Bt est la transposée de B

    // Résoudre By = Atb pour y par substitution avant
    vector<double> y = substitution_avant(B, Atb);
    // Résoudre B^Tx = y pour x par substitution arrière
    vector<double> x_chapeau = substitution_arriere(Bt, y);

    return x_chapeau;
}

//Test :
// Définissons une matrice A et un vecteur b
    // vector<vector<double>> A = {
    //     {2, -1, 0, 0},
    //     {-1, 2, -1, 0},
    //     {0, -1, 2, -1},
    //     {0, 0, -1, 2}
    // };
    // vector<double> b = {1, 0, 0, 0};

    // // Utilisons la fonction pour calculer la solution des moindres carrés
    // try {
    //     vector<double> x_chapeau = resoudre_equation_normale(A, b);
        
    //     // Affichons le résultat
    //     cout << "La solution x_chapeau est:" << endl;
    //     for (double val : x_chapeau) {
    //         cout << val << " ";
    //     }
    //     cout << endl;
    // } catch (const std::exception& e) {
    //     cerr << "Erreur: " << e.what() << endl;
    // }

    // return 0;

//******************************************************************************************************************************

void readHousingData(const string& filename, vector<vector<double>>& A, vector<double>& b,
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

        vector<double> rowA{1.0};  // Start with the intercept term

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
    // vector<vector<double>> A;
    // vector<double> b;
    // string filename = "housing.data.txt";

    // readHousingData(filename, A, b);

    // vector<double> x_chapeau = resoudre_equation_normale(A, b);

    //     // Affichons le résultat
    //     cout << "La solution x_chapeau est:" << endl;
    //     for (double val : x_chapeau) {
    //         cout << val << " ";
    //     }
    //     cout << endl;


    // return 0;

//******************************************************************************************************************************

double calculeSigmaChapeau(const vector<vector<double>>& A, const vector<double>& b, const vector<double>& theta) {
    vector<double> residus(b.size());
    double sommeDesCarresDesResidus = 0.0;

    // Calculer les valeurs prédites
    for (size_t i = 0; i < A.size(); ++i) {
        double predite = 0.0;
        for (size_t j = 0; j < A[i].size(); ++j) {
            predite += A[i][j] * theta[j];
        }
        residus[i] = b[i] - predite;
        sommeDesCarresDesResidus += std::pow(residus[i], 2);
    }

    // Calculer sigma chapeau
    return std::sqrt(sommeDesCarresDesResidus / b.size());
}

    // vector<vector<double>> A;
    // vector<double> b;
    // string filename = "housing.data.txt"; // Ensure this path is correct

    // readHousingData(filename, A, b);
    // vector<double> theta = resoudre_equation_normale(A, b);
    // double sigmaHat = calculeSigmaChapeau(A, b, theta);

    // cout << "Sigma Chapeau: " << sigmaHat << endl;

    // return 0; 

//******************************************************************************************************************************
// For debug purposes
void writeMatrixToFile(const vector<vector<double>>& matrix, const string& filename) {
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
void generatePolynomialTerms(const vector<double>& vars, int degree, vector<double>& row) {
    row.push_back(1.0);  // Add intercept term
    for (int i = 1; i <= degree; ++i) {
        for (double var : vars) {
            row.push_back(pow(var, i));  // Add variable raised to the power 'i'
        }
    }
}

// Main function to read data and create the matrix A and vector b
void createPolynomialRegressionMatrix(vector<vector<double>>& data, int degree, vector<vector<double>>& A, vector<double>& b) {
    for (const auto& row : data) {
        vector<double> A_row;
        vector<double> vars(row.begin(), row.end() - 1);  // Extract all variables (assuming last column is the dependent variable)
        generatePolynomialTerms(vars, degree, A_row);
        A.push_back(A_row);
        b.push_back(row.back());  // Last element is the dependent variable
    }
}






int main(){

    // Example usage:
    // Suppose data is a vector of vectors, each containing [x1, x2, ..., xn, y]
    // For example, data for a model with 2 predictors might look like this:
    vector<vector<double>> data = {
        {1.0, 2.0, 10.0},
        {1.5, 3.5, 15.0},
        {2.0, 4.0, 20.0}
    };

    vector<vector<double>> A;
    vector<double> b;
    int polynomialDegree = 4;  // Polynomial degree

    createPolynomialRegressionMatrix(data, polynomialDegree, A, b);

    // Print matrix A and vector b
    cout << "Matrix A:" << endl;
    for (const auto& row : A) {
        for (auto elem : row) cout << elem << " ";
        cout << endl;
    }

    cout << "Vector b:" << endl;
    for (auto elem : b) cout << elem << " ";
    cout << endl;

    return 0;

}

    // vector<vector<double>> A;
    // vector<double> b;
    // readHousingData("housing.data.txt", A, b, true, true, true, false, 0);  
    // vector<double> theta = resoudre_equation_normale(A, b);
    // double sigmaHat = calculeSigmaChapeau(A, b, theta);

    // cout << "Sigma Chapeau: " << sigmaHat << endl;

    // writeMatrixToFile(A, "outputMatrix.txt");  // Specify the output filename

    // return 0; 