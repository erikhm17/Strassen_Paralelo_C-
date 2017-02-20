#include "mpi.h"
#include <string>
#include <stdio.h>
#include <iostream>
#include <math.h> 


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
double ** sub(double **A, double ** B, int size)
{
	int n = size;
	double ** C = newMatriz(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			C[i][j] = A[i][j] - B[i][j];
	return C;
}
/** Funcion para sumar dos matrices  **/
double ** add(double ** A, double** B, int size)
{
	int n = size;
	double ** C = newMatriz(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			C[i][j] = A[i][j] + B[i][j];
	return C;
}
/** Funcion para separar la matriz padre dentro de los hijos**/
void fraccionar(double** P, double** C, int iB, int jB, int size)
{
	for (int i1 = 0, i2 = iB; i1 < size; i1++, i2++)
		for (int j1 = 0, j2 = jB; j1 < size; j1++, j2++)
			C[i1][j1] = P[i2][j2];
}

/** Funcion para unir las matrices hijas y matriz padre **/
void juntar(double** C, double** P, int iB, int jB, int size)
{
	for (int i1 = 0, i2 = iB; i1 < size; i1++, i2++)
		for (int j1 = 0, j2 = jB; j1 < size; j1++, j2++)
			P[i2][j2] = C[i1][j1];
}

void escribirArreglo(double * arreglo, int longitud) {
	for (int i = 0; i < longitud; i++){
		cout << " "<<arreglo[i];
	}
	cout << "\n";
}

void escribirMatriz(double **matriz, int size) {
	for (int k = 0; k < size; k++) {
		for (int i = 0; i < size; i++) {
			cout << " " << matriz[k][i];
			//printf("%f", matriz[k][i]);
		}
		cout << "\n ";
		//printf("\n");
	}
}
double * convertirEnArreglo(double **matriz, int dimensionMatriz) {
	int sizeArray= dimensionMatriz*dimensionMatriz; // para que sea lineal es N*N 
	double *arreglo = new double[sizeArray];
	int count = 0;
	for (int i = 0; i < dimensionMatriz; i++){
		for (int j  = 0; j < dimensionMatriz; j++)	{
			arreglo[count] = matriz[i][j];
			count++;
		}
	}
	return arreglo;
}

double ** convertirEnMatriz(double *arreglo,int dimensionArreglo) {
	int sizeMatriz = sqrt(dimensionArreglo); /* calculamos la dimension de la matriz*/
	double **matriz = newMatriz(sizeMatriz,sizeMatriz);
	int count = 0;
	for (int i = 0; i < sizeMatriz; i++){
		for (int j = 0; j < sizeMatriz; j++) {
			matriz[i][j] = arreglo[count];
			count++;
		}
	}
	return matriz;
}
double **sumarMatrices(double ** matrizA , double **matrizB,int dimension) {
	double ** matrizSuma = newMatriz(dimension,dimension);
	for (int i = 0; i < dimension; i++) {
		for (int j = 0; j < dimension; j++) {
			matrizSuma[i][j] = matrizA[i][j]+matrizB[i][j];
		}
	}
	return matrizSuma;
}
double **restarMatrices(double ** matrizA, double **matrizB, int dimension) {
	double ** matrizResta = newMatriz(dimension, dimension);
	for (int i = 0; i < dimension; i++) {
		for (int j = 0; j < dimension; j++) {
			matrizResta[i][j] = matrizA[i][j] - matrizB[i][j];
		}
	}
	return matrizResta;
}

class Strassen {
public:
	Strassen() {

	}
	/** funcion para sumar matrices */
	double **multiply(double **A, double **B, int size, int *rank, int *sizeProcesador) {
		int n = size;
		double ** R = newMatriz(n, n);
		/* caso base */
		if (n == 1) {
			R[0][0] = A[0][0] * B[0][0];
		}
		else {
			double ** A11 = newMatriz(n / 2, n / 2);
			double ** A12 = newMatriz(n / 2, n / 2);
			double ** A21 = newMatriz(n / 2, n / 2);
			double ** A22 = newMatriz(n / 2, n / 2);
			double ** B11 = newMatriz(n / 2, n / 2);
			double ** B12 = newMatriz(n / 2, n / 2);
			double ** B21 = newMatriz(n / 2, n / 2);
			double ** B22 = newMatriz(n / 2, n / 2);

			/** Dividing matrix A into 4 halves **/
			/** Dividiendo la matriz A en 4 particiones, (infinitamente ... ) */
			fraccionar(A, A11, 0, 0, n / 2);
			fraccionar(A, A12, 0, n / 2, n / 2);
			fraccionar(A, A21, n / 2, 0, n / 2);
			fraccionar(A, A22, n / 2, n / 2, n / 2);
			/** Dividing matrix B into 4 halves **/
			fraccionar(B, B11, 0, 0, n / 2);
			fraccionar(B, B12, 0, n / 2, n / 2);
			fraccionar(B, B21, n / 2, 0, n / 2);
			fraccionar(B, B22, n / 2, n / 2, n / 2);

			/**
			M1 = (A11 + A22)(B11 + B22)
			M2 = (A21 + A22) B11
			M3 = A11 (B12 - B22)
			M4 = A22 (B21 - B11)
			M5 = (A11 + A12) B22
			M6 = (A21 - A11) (B11 + B12)
			M7 = (A12 - A22) (B21 + B22)
			**/

			double ** M1 = multiply(add(A11, A22, n / 2), add(B11, B22, n / 2), n / 2, rank, sizeProcesador);
			double ** M2 = multiply(add(A21, A22, n / 2), B11, n / 2, rank, sizeProcesador);
			double ** M3 = multiply(A11, sub(B12, B22, n / 2), n / 2, rank, sizeProcesador);
			double ** M4 = multiply(A22, sub(B21, B11, n / 2), n / 2, rank, sizeProcesador);
			double ** M5 = multiply(add(A11, A12, n / 2), B22, n / 2, rank, sizeProcesador);
			double ** M6 = multiply(sub(A21, A11, n / 2), add(B11, B12, n / 2), n / 2, rank, sizeProcesador);
			double ** M7 = multiply(sub(A12, A22, n / 2), add(B21, B22, n / 2), n / 2, rank, sizeProcesador);

			/**
			C11 = M1 + M4 - M5 + M7
			C12 = M3 + M5
			C21 = M2 + M4
			C22 = M1 - M2 + M3 + M6
			**/
			double ** C11 = add(sub(add(M1, M4, n / 2), M5, n / 2), M7, n / 2);
			double ** C12 = add(M3, M5, n / 2);
			double ** C21 = add(M2, M4, n / 2);
			double ** C22 = add(sub(add(M1, M3, n / 2), M2, n / 2), M6, n / 2);

			/** join 4 halves into one result matrix **/
			juntar(C11, R, 0, 0, n / 2);
			juntar(C12, R, 0, n / 2, n / 2);
			juntar(C21, R, n / 2, 0, n / 2);
			juntar(C22, R, n / 2, n / 2, n / 2);
		}

		return  R;


	}
};
int main(int argc, char **argv) {
	int rank, size;
	MPI_Status status;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // los ides de cada proceso 
	MPI_Comm_size(MPI_COMM_WORLD, &size); // numero de procesos que lo componen
	
										  /* ponemos una dimension de 4*/
	int dimension = 32;
	
	int mitad = dimension/2;
	double ** matrizA = newMatriz(dimension, dimension);
	double ** matrizB = newMatriz(dimension, dimension);
	

	if (rank==0){

		for (int k = 0; k < dimension; k++) {
			for (int u = 0; u < dimension; u++) {
				matrizA[k][u] =1; // generamos numeros
			}
		}
		for (int k = 0; k < dimension; k++) {
			for (int u = 0; u < dimension; u++) {
				matrizB[k][u] = 1; // generamos numeros
			}
		}

		cout << "Matriz A \n ";
		escribirMatriz(matrizA, dimension);
		cout << "\n ";
		cout << "Matriz B \n ";
		escribirMatriz(matrizB, dimension);
		cout << "\n ";
	}
		
		/*enviar todo los As */
		double *arrayMatrizA = convertirEnArreglo(matrizA,dimension);
		double *arrayMatrizB = convertirEnArreglo(matrizB, dimension);

		MPI_Bcast(arrayMatrizA,dimension*dimension, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Bcast(arrayMatrizB, dimension*dimension, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		//printf("proceso %d, tiene la siguiente matriz :\n", rank);
		matrizA = convertirEnMatriz(arrayMatrizA, dimension*dimension);
		//escribirMatriz(matrizA, dimension);
		//cout << "\n";

		matrizB = convertirEnMatriz(arrayMatrizB, dimension*dimension);
		//escribirMatriz(matrizB, dimension);
		Strassen *s = new Strassen();
		int n = dimension/2; // n es la sub_dimension 
		double *M1_Array = new double[n*n];
		if (rank == 0) {
			double ** A11 = newMatriz(n, n);
			double ** A22 = newMatriz(n, n);
			double ** B11 = newMatriz(n, n);
			double ** B22 = newMatriz(n, n);

			fraccionar(matrizA, A11, 0, 0, n);
			fraccionar(matrizA, A22, n, n, n);
			fraccionar(matrizB, B11, 0, 0, n);
			fraccionar(matrizB, B22, n, n, n);

			double ** T1 = sumarMatrices(A11, A22, n);
			double ** T2 = sumarMatrices(B11, B22, n);
			double **M1_Aux = s->multiply(T1, T2, n, &rank, &size); // M1
			M1_Array = convertirEnArreglo(M1_Aux, n); // M1_Aux_Envio seria de n x n
	/*		cout << " ##############################  : \n ";

			cout << " T1 : \n ";
			escribirMatriz(T1, n);
			cout << " \n ";

			cout << " T2 : \n ";
			escribirMatriz(T2, n);
			cout << " \n ";

			cout << "  B11  : \n ";
			escribirMatriz(B11, n);
			cout << "  B22 \n ";
			escribirMatriz(B22, n);
			cout << "  SUMA es T2:  \n ";
			*/

			/*
			
			cout << " M1 : \n ";
			escribirMatriz(M1_Aux,n);
			cout << " \n ";
			
			
			*/

			/* MPI_Send : no puede enviarse a si mismo */
			//MPI_Send(M1_Aux_Envio, n*n, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD) // mal;

		}
		if (rank == 1) {
			double ** A11 = newMatriz(n, n);
			double ** A21 = newMatriz(n, n);
			double ** A22 = newMatriz(n, n);
			double ** B11 = newMatriz(n, n);
			double ** B12 = newMatriz(n, n);
			double ** B22 = newMatriz(n, n);

			fraccionar(matrizA, A11, 0, 0, n);
			fraccionar(matrizA, A21, n, 0, n);
			fraccionar(matrizA, A22, n, n, n);
			fraccionar(matrizB, B11, 0, 0, n);
			fraccionar(matrizB, B12, 0, n, n);
			fraccionar(matrizB, B22, n, n, n);

			double ** T3 = sumarMatrices(A21, A22, n);
			double ** M2_Aux = s->multiply(T3, B11, n, &rank, &size); // M2
			double * M2_Aux_Envio = convertirEnArreglo(M2_Aux, n);

			double **T4 = restarMatrices(B12, B22, n);
			double **M3_Aux = s->multiply(A11, T4, n, &rank, &size); // M3
			double * M3_Aux_Envio = convertirEnArreglo(M3_Aux, n);

/*			cout << " ##############################  : \n ";

			cout << " T3 : \n ";
			escribirMatriz(T3, n);
			cout << " \n ";

			cout << " T4 : \n ";
			escribirMatriz(T4, n);
			cout << " \n ";
			*/

			/*


			cout << " M2 : \n ";
			escribirMatriz(M2_Aux, n);
			cout << " M3 \n ";
			escribirMatriz(M3_Aux, n);
			cout << " \n ";


			*/

			//double * mensaje = new double[1];
			//mensaje[0] = 666;
			
			
			/*
			
			cout << "Arreglo a enviar : ";
			escribirArreglo(M2_Aux_Envio,n*n);
			cout << "\n ";


			*/


			MPI_Send(M2_Aux_Envio, n*n, MPI_DOUBLE, 0, 200, MPI_COMM_WORLD);
			MPI_Send(M3_Aux_Envio, n*n, MPI_DOUBLE, 0, 300, MPI_COMM_WORLD);

		}
		if (rank == 2) {
	
			double ** A11 = newMatriz(n, n);
			double ** A12 = newMatriz(n, n);
			double ** A21 = newMatriz(n, n);
			double ** A22 = newMatriz(n, n);
			double ** B11 = newMatriz(n, n);
			double ** B12 = newMatriz(n, n);
			double ** B21 = newMatriz(n, n);
			double ** B22 = newMatriz(n, n);

			fraccionar(matrizA, A11, 0, 0, n);
			fraccionar(matrizA, A12, 0, n, n);
			fraccionar(matrizA, A22, n, n, n);
			/** Dividing matrix B into 4 halves **/
			fraccionar(matrizB, B11, 0, 0, n);
			fraccionar(matrizB, B21, n, 0, n);
			fraccionar(matrizB, B22, n, n, n);
			double ** T5 = restarMatrices(B21, B11, n);
			double ** T6 = sumarMatrices(A11, A12, n);
			
			double ** M4_Aux = s->multiply(A22, T5, n, &rank, &size); // M4
			double * M4_Aux_Envio = convertirEnArreglo(M4_Aux, n);
			
			double ** M5_Aux = s->multiply(T6, B22, n, &rank, &size); // M5
			double * M5_Aux_Envio = convertirEnArreglo(M5_Aux, n);

			MPI_Send(M4_Aux_Envio, n*n, MPI_DOUBLE, 0, 400, MPI_COMM_WORLD);
			MPI_Send(M5_Aux_Envio, n*n, MPI_DOUBLE, 0, 500, MPI_COMM_WORLD);

		/*	cout << " ##############################  : \n ";

			cout << " T5 : \n ";
			escribirMatriz(T5, n);
			cout << " \n ";

			cout << " T6 : \n ";
			escribirMatriz(T6, n);
			cout << " \n ";
			*/
			/*

			cout << " M4 : \n ";
			escribirMatriz(M4_Aux, n);
			cout << " M5 :  \n ";
			escribirMatriz(M5_Aux, n);
			cout << " \n ";


			*/

		}
		if (rank == 3) {
			
			double ** A11 = newMatriz(n, n);
			double ** A12 = newMatriz(n, n);
			double ** A21 = newMatriz(n, n);
			double ** A22 = newMatriz(n, n);
			double ** B11 = newMatriz(n, n);
			double ** B12 = newMatriz(n, n);
			double ** B21 = newMatriz(n, n);
			double ** B22 = newMatriz(n, n);

			fraccionar(matrizA, A11, 0, 0, n);
			fraccionar(matrizA, A12, 0, n, n);
			fraccionar(matrizA, A21, n, 0, n);
			fraccionar(matrizA, A22, n, n, n);
			/** Dividing matrix B into 4 halves **/
			fraccionar(matrizB, B11, 0, 0, n);
			fraccionar(matrizB, B12, 0, n, n);
			fraccionar(matrizB, B21, n, 0, n);
			fraccionar(matrizB, B22, n, n, n);

			double ** T7 = restarMatrices(A21, A11, n);
			double ** T8 = sumarMatrices(B11, B12, n);
			double ** M6_Aux = s->multiply(T7, T8, n, &rank, &size); // M6

			double ** T9 = restarMatrices(A12, A22, n);
			double ** T10 = sumarMatrices(B21, B22, n);
			double ** M7_Aux = s->multiply(T9, T10, n, &rank, &size); // M7


			double * M6_Aux_Envio = convertirEnArreglo(M6_Aux, n);
			double * M7_Aux_Envio = convertirEnArreglo(M7_Aux, n);

			MPI_Send(M6_Aux_Envio, n*n, MPI_DOUBLE, 0, 600, MPI_COMM_WORLD);
			MPI_Send(M7_Aux_Envio, n*n, MPI_DOUBLE, 0, 700, MPI_COMM_WORLD);

			cout << " ##############################  : \n ";

		/*	cout << " T7 : \n ";
			escribirMatriz(T7, n);
			cout << " \n ";

			cout << " T8 : \n ";
			escribirMatriz(T8, n);
			cout << " \n ";

			cout << " T9 : \n ";
			escribirMatriz(T9, n);
			cout << " \n ";

			cout << " T10 : \n ";
			escribirMatriz(T10, n);
			cout << " \n ";
			*/

			/*

			cout << " M6 : \n ";
			escribirMatriz(M6_Aux, n);
			cout << " M7 :\n ";
			escribirMatriz(M7_Aux, n);
			cout << " \n ";

			*/

		}
		if (rank==0){
			/*creo mi arreglo de Ms*/
			
			double *M2_Array = new double[n*n];
			double *M3_Array = new double[n*n];

			double *M4_Array = new double[n*n];
			double *M5_Array = new double[n*n];
			
			double *M6_Array = new double[n*n];
			double *M7_Array = new double[n*n];
			

			/** el M1_Array ya se le asigna localmente no es necesario un envio */
			
			MPI_Recv(M2_Array, n*n, MPI_DOUBLE, 1, 200, MPI_COMM_WORLD, &status);
			MPI_Recv(M3_Array, n*n, MPI_DOUBLE, 1, 300, MPI_COMM_WORLD, &status);
			MPI_Recv(M4_Array, n*n, MPI_DOUBLE, 2, 400, MPI_COMM_WORLD, &status);
			MPI_Recv(M5_Array, n*n, MPI_DOUBLE, 2, 500, MPI_COMM_WORLD, &status);
			MPI_Recv(M6_Array, n*n, MPI_DOUBLE, 3, 600, MPI_COMM_WORLD, &status);
			MPI_Recv(M7_Array, n*n, MPI_DOUBLE, 3, 700, MPI_COMM_WORLD, &status);

			//double ** M1_Aux_matriz = convertirEnMatriz(M1, n*n);

			cout << "M1 es : ";
			escribirArreglo(M1_Array, n*n);
			cout << "M2 es : ";
			escribirArreglo(M2_Array, n*n);
			cout << "M3 es : ";
			escribirArreglo(M3_Array, n*n);
			cout << "M4 es : ";
			escribirArreglo(M4_Array, n*n);
			cout << "M5 es : ";
			escribirArreglo(M5_Array, n*n);
			cout << "M6 es : ";
			escribirArreglo(M6_Array, n*n);
			cout << "M7 es : ";
			escribirArreglo(M7_Array, n*n);
	
			double ** M1 = convertirEnMatriz(M1_Array, n*n);
			double ** M2 = convertirEnMatriz(M2_Array, n*n);
			double ** M3 = convertirEnMatriz(M3_Array, n*n);
			double ** M4 = convertirEnMatriz(M4_Array, n*n);
			double ** M5 = convertirEnMatriz(M5_Array, n*n);
			double ** M6 = convertirEnMatriz(M6_Array, n*n);
			double ** M7 = convertirEnMatriz(M7_Array, n*n);

			double ** C11 = add(sub(add(M1, M4, n), M5, n), M7, n);
			double ** C12 = add(M3, M5, n);
			double ** C21 = add(M2, M4, n);
			double ** C22 = add(sub(add(M1, M3, n), M2, n ), M6, n );
			
			double **R = newMatriz(dimension,dimension);
			/** join 4 halves into one result matrix **/
			juntar(C11, R, 0, 0, n);
			juntar(C12, R, 0, n , n );
			juntar(C21, R, n, 0, n );
			juntar(C22, R, n , n , n );

			cout << " La matriz resultante es : \n";
			escribirMatriz(R,dimension);
		}
		MPI_Finalize();
		//system("PAUSE");
		//cin.get();
	

	
	return 0;
}