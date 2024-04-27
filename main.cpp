#include <iostream>
#include <string>
#include <vector>   
#include <cmath> 

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


bool isValide(vector<vector<double>> matrice){
    for (vector<double> e: matrice){
        if(matrice[0].size() != e.size()){
            return 0;
        }
    }return 1;
}




//******************************************************************************************************************************
int main(){

    vector<vector<double>> matrix = {
        {1, 2, 3},
        {4, 5, 6}
    };

    vector<vector<double>> transposed = transpose_matrice(matrix);

    for (const auto& row : transposed) {
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }

    return 0;
}