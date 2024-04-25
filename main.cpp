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
int main(){
    vector<double> v1,v2;
    v1.push_back(1.0);
    v1.push_back(6.0);

    v2.push_back(1.0);
    v2.push_back(9.0);

    vector<double> result = addition(v1,v2);
    for(double val : result) {
        cout << val << " ";
    }
}