
#include <string>
#include <stdio.h>
#include <iostream>

using namespace std;

double **newMatriz(int filas, int columnas) {
    double **matriz = new double*[filas];
    for (int i = 0; i < filas; i++)
    {
        matriz[i] = new double[columnas];
    }
    return matriz;
}
/** Funtion to sub two matrices **/
double ** sub(double **A, double ** B,int size)
{
    int n = size;
    double ** C = newMatriz(n,n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = A[i][j] - B[i][j];
    return C;
}
/** Funcion para sumar dos matrices  **/
double ** add(double ** A, double** B,int size)
{
    int n = size;
    double ** C = newMatriz(n,n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = A[i][j] + B[i][j];
    return C;
}
/** Funcion para separar la matriz padre dentro de los hijos**/
void split(double** P, double** C, int iB, int jB,int size)
{
    for(int i1 = 0, i2 = iB; i1 < size; i1++, i2++)
        for(int j1 = 0, j2 = jB; j1 < size; j1++, j2++)
            C[i1][j1] = P[i2][j2];
}

/** Funcion para unir las matrices hijas y matriz padre **/
void join(double** C, double** P, int iB, int jB,int size)
{
    for(int i1 = 0, i2 = iB; i1 < size; i1++, i2++)
        for(int j1 = 0, j2 = jB; j1 < size; j1++, j2++)
            P[i2][j2] = C[i1][j1];
}
void escribirMatriz(double **matriz, int size){
    for (int k = 0; k < size; k++) {
        for (int i = 0; i < size; i++) {
            cout<<" "<<matriz[k][i];
        }
        cout<<"\n ";
    }
}
class Strassen{
public:
    Strassen(){

    }
    /** funcion para sumar matrices */
    double **multiply(double **A, double **B, int size){
        int n = size;
        double ** R= newMatriz(n, n);
        /* caso base */
        if (n == 1){
            R[0][0] = A[0][0] * B[0][0];
        }else{
            double ** A11 = newMatriz(n/2,n/2);
            double ** A12 = newMatriz(n/2,n/2);
            double ** A21 = newMatriz(n/2,n/2);
            double ** A22 = newMatriz(n/2,n/2);
            double ** B11 = newMatriz(n/2,n/2);
            double ** B12 = newMatriz(n/2,n/2);
            double ** B21 = newMatriz(n/2,n/2);
            double ** B22 = newMatriz(n/2,n/2);

            /** Dividing matrix A into 4 halves **/
            /** Dividiendo la matriz A en 4 particiones, (infinitamente ... ) */
            split(A, A11, 0 , 0,n/2);
            split(A, A12, 0 , n/2,n/2);
            split(A, A21, n/2, 0,n/2);
            split(A, A22, n/2, n/2,n/2);
            /** Dividing matrix B into 4 halves **/
            split(B, B11, 0 , 0,n/2);
            split(B, B12, 0 , n/2,n/2);
            split(B, B21, n/2, 0,n/2);
            split(B, B22, n/2, n/2,n/2);

            /**
             M1 = (A11 + A22)(B11 + B22)
             M2 = (A21 + A22) B11
             M3 = A11 (B12 - B22)
             M4 = A22 (B21 - B11)
             M5 = (A11 + A12) B22
             M6 = (A21 - A11) (B11 + B12)
             M7 = (A12 - A22) (B21 + B22)
             **/

            double ** M1 = multiply(add(A11, A22,n/2), add(B11, B22,n/2),n/2);
            double ** M2 = multiply(add(A21, A22,n/2), B11,n/2);
            double ** M3 = multiply(A11, sub(B12, B22,n/2),n/2);
            double ** M4 = multiply(A22, sub(B21, B11,n/2),n/2);
            double ** M5 = multiply(add(A11, A12,n/2), B22,n/2);
            double ** M6 = multiply(sub(A21, A11,n/2), add(B11, B12,n/2),n/2);
            double ** M7 = multiply(sub(A12, A22,n/2), add(B21, B22,n/2),n/2);

            /**
             C11 = M1 + M4 - M5 + M7
             C12 = M3 + M5
             C21 = M2 + M4
             C22 = M1 - M2 + M3 + M6
             **/
            double ** C11 = add(sub(add(M1, M4,n/2), M5,n/2), M7,n/2);
            double ** C12 = add(M3, M5,n/2);
            double ** C21 = add(M2, M4,n/2);
            double ** C22 = add(sub(add(M1, M3,n/2), M2,n/2), M6,n/2);

            /** join 4 halves into one result matrix **/
            join(C11, R, 0 , 0,n/2);
            join(C12, R, 0 , n/2,n/2);
            join(C21, R, n/2, 0,n/2);
            join(C22, R, n/2, n/2,n/2);
        }

        return  R;


    }
};

class persona{
    int dni;
    string nombre ;
    string apellido;
    int edad;

    public :persona(int dni,string nombre, string apellido ,int edad){
        this->dni = dni;
        this->nombre = nombre;
        this->apellido = apellido;
        this->edad = edad;

    }
    void escribir() {
        cout<<"Nombre: "<<this->nombre<<" "
            " apellido : "<<this->apellido<<" edad : "<<this->edad<<"\n";
    }
    /*
     * persona *p = new persona(123123,"Erik","choqque",24);
       p->escribir();
     * */
};

int main(){

    printf(" Strassen Multiplication Algorithm Test \n");
    /** Make an object of Strassen class **/
    int N ;

    Strassen *s = new Strassen();
    cout<<"Ingrese el orden N \n";
    cin>>N;

    cout<<"Ingrese valores de la matriz A \n";
    double **A = newMatriz(N,N);
    int auxA,valoresA;

    for (int i = 0; i<N;i++){
        for (int j = 0; j < N; j++) {
            cin>>valoresA;
            auxA= valoresA;
            A[i][j] = auxA;
        }
    }

    cout<<"Ingrese valores de la matriz B \n";
    double **B = newMatriz(N,N);
    int auxB,valoresB;

    for (int i = 0; i<N;i++){
        for (int j = 0; j < N; j++) {
            cin>>valoresB;
            auxB= valoresB;
            B[i][j] = auxB;
        }
    }
    escribirMatriz(A,N);
    cout<<"\n ";
    escribirMatriz(B,N);

    double **C = s->multiply(A,B,N);
    cout<<" Producto de matrices entre A y B es : \n";
    escribirMatriz(C,N);

    return 0;
}