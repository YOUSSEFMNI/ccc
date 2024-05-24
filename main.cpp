#include <iostream>
#include "Matrix.h"

void executeMethod(const Matrix& matrix, const std::string& method) {
    try {
        Matrix inverseMatrix = (method == "gauss") ? matrix.inverseGauss() : matrix.inverseCramer();
        std::cout << "Inverse de la matrice en utilisant : " << ((method == "gauss") ? "Gauss" : "Cramer") << ":" << std::endl;
        inverseMatrix.display();

        // verifier le resultat
        Matrix product = matrix.multiply(inverseMatrix);
        std::cout << "Produit de la matrice originale et son inverse :" << std::endl;
        product.roundSmallValues();
        product.display();

        if (product.isIdentity()) {
            std::cout << "L'invertion est correcte , le produit est une matrice identitee " << std::endl;
        }
        else {
            std::cout << "l'invertion n'est pas correcte le produit n'est pas une matrice identitee " << std::endl;
        }
    }
    catch (const std::logic_error& e) {
        std::cerr << "Erreur : " << e.what() << std::endl;
    }
}

int main() {
    int rows, cols;
    std::cout << "Entrer le nombre de lignes ";
    std::cin >> rows;
    std::cout << "Entrer le nombre de collones ";
    std::cin >> cols;

    Matrix matrix(rows, cols);
    matrix.input();
    std::cout << "matrice originale:" << std::endl;
    matrix.display();

    if (!matrix.isSquare()) {
        std::cout << "Ce n'est pas une matrice carrée." << std::endl;
        char choix;
        std::cout << "Voulez-vous effectuer la pseudo-inversion de la matrice ? (o/n): ";
        std::cin >> choix;
        if (choix == 'o' || choix == 'O') {
            Matrix A_pseudoInv = matrix.pseudoInverse();
            if (A_pseudoInv.getRows() != 0 && A_pseudoInv.getCols() != 0) {
                std::cout << "Pseudo-inverse de la matrice A :" << std::endl;
                A_pseudoInv.display();
                Matrix product = matrix.multiply(A_pseudoInv);
                std::cout << "Produit de la matrice originale et son inverse :" << std::endl;
                product.roundSmallValues();
                product.display();
            }
            else {
                std::cout << "La pseudo-inverse n'a pas pu être calculée." << std::endl;
            }
        }
    }
    else {
        std::cout << "C'est une matrice carrée. Vous pouvez effectuer l'inversion normale." << std::endl;
        std::cout << "Coisir votre methode (gauss/cramer): ";
        std::string method;
        std::cin >> method;

        executeMethod(matrix, method);
    }

    return 0;
}