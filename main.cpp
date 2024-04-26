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

vector<vector<double>> transpose_matrice(const vector<vector<double>>& matrix) {
    if (matrix.empty()) return {};

    // Create a new matrix with the dimensions swapped
    vector<vector<double>> transposed(matrix[0].size(), vector<double>(matrix.size()));

    // Fill the transposed matrix
    for (size_t i = 0; i < matrix.size(); i++) {
        for (size_t j = 0; j < matrix[0].size(); j++) {
            transposed[j][i] = matrix[i][j];
        }
    }

    return transposed;
}
//******************************************************************************************************************************

vector<double> produit_matrices(){
    
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