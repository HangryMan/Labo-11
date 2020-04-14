#include <stdio.h>
#include <math.h>

#define M 3         //размерность входящего массива
#define e 0.0000001     //точность вычисления

void convertToIterForm(void);
void printMatrix(void);
void convertToDiagPrev(void);
int convergenceCond(void);
void iterProcessingJacobi(void);
void iterProcessingZeidel(void);
void printMatrixX(void);
void printRes(void);

float matrix[M][M+1] = {  {0.12, -0.43, 0.14, -0.17},
                          {-0.07, 0.34, 0.72, 0.62},
                          {1.18, -0.08, -0.25, 1.12}};
float matrixX[M] = {0}, bufMatrixX[M] = {0};

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<MAIN
int main(void) {
    puts("\tПервичная матрица:");
    printMatrix();

    puts(" Приведённая к диаг преобладанию:");
    convertToDiagPrev();

    puts("    Приведённая к итер форме:");
    convertToIterForm();

    puts("\tПроверка на сходимость:");
    if (convergenceCond())
        puts("Success: convergenceCond return 1\n");
    else {
        puts("Error: convergenceCond return 0\n");
        fflush(stdin);
        getchar();
        return 0;
    }

    printf("Выбор метода решения:\n1)Якоби\n2)Зейделя\n>");
    if (getchar() == '1') iterProcessingJacobi();
    else iterProcessingZeidel();

    printRes();
    return 0;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void printMatrix(void) {
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M+1; ++j) {
            printf("%8.4f", matrix[i][j]);
        }
        puts("");
    }
    puts("");
}
void convertToDiagPrev(void) {
    float buff = 0;
    unsigned short indexes[M] = {0};


    for (int i = 0; i < M; ++i) {
        float maxAbsolValue = fabsf(matrix[i][0]);

        for (int j = 0; j < M; ++j) {
            if (fabsf(matrix[i][j]) > maxAbsolValue) {
                maxAbsolValue = fabsf(matrix[i][j]);
                indexes[i] = j;
            }
        }
    }

    for (int i = 0; i < M; ++i) {
        float min = indexes[0];
        unsigned short minInx = 0;

        for (int k = i; k < M; ++k) {
            if (indexes[k] < min) {
                min = indexes[k];
                minInx = k;
            }
        }

        for (int j = 0; j < M+1; ++j) {
            buff = matrix[i][j];
            matrix[i][j] = matrix[indexes[minInx]][j];
            matrix[indexes[minInx]][j] = buff;
        }

        printMatrix();
    }
}
void convertToIterForm(void) {
    for (int i = 0; i < M; ++i) {
        float diagElem = matrix[i][i];
        matrix[i][i] = 0;
        for (int j = 0; j < M+1; ++j) {
            matrix[i][j] = (j < M? -1: 1)*matrix[i][j]/diagElem;
        }
        printMatrix();
    }
}
int convergenceCond(void) {
    float sum = 0;
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            sum += pow(matrix[i][j], 2);
        }
    }
    //printf("%f\n", sqrt(sum));
    if (sqrt(sum) < 1) return 1; else return 0;
}
void iterProcessingJacobi(void) {
    puts("\n\tИтеррации:");

    for (int i = 0; i < M; ++i)         //Нулевая итеррация (X(0) = Beta)
        bufMatrixX[i] = matrix[i][M];

    unsigned int counter = 0;
    while (counter < M) {
        counter = 0;
        for (int i = 0; i < M; ++i) {
            float delta = 0;
            for (int j = 0; j < M; ++j)
                delta += matrix[i][j] * bufMatrixX[j]; // delta - это alfa*X(k-1)

            matrixX[i] = matrix[i][M] + delta;
            if (fabsf(matrixX[i] - bufMatrixX[i]) <= e) {
                counter++;
                continue;}
        }

        for (int i = 0; i < M; ++i) {
            bufMatrixX[i] = matrixX[i];
        }

        printMatrixX();
    }
}
void iterProcessingZeidel(void) {
    puts("\n\tИтеррации:");

    for (int i = 0; i < M; ++i)         //Нулевая итеррация (X(0) = Beta)
        bufMatrixX[i] = matrix[i][M];

    unsigned int counter = 0;
    while (counter < M) {
        counter = 0;
        for (int i = 0; i < M; ++i) {
            float delta1 = 0;
            for (int j = i - 1; j < M; ++j)
                delta1 += matrix[i][j] * bufMatrixX[j]; // delta1 - это alfa*X(k-1)

            matrixX[i] = matrix[i][M] + delta1;

            float delta2 = 0;
            for (int j = 0; j < i - 1; ++j)
                delta2 += matrix[i][j] * matrixX[j]; // delta2 - это alfa*X(k)

            matrixX[i] += delta2;

            if (fabsf(matrixX[i] - bufMatrixX[i]) <= e) {
                counter++;
                continue;}
        }

        for (int i = 0; i < M; ++i) {
            bufMatrixX[i] = matrixX[i];
        }

        printMatrixX();
    }
}
void printMatrixX(void) {
    static num = 1;
    printf("[%d] ", num++);
    for (int i = 0; i < M; ++i) {
        printf("%f ", matrixX[i]);
    }
    puts("");
}
void printRes(void) {
    for (int i = 0; i <M; ++i) {
        printf("X%d = %f\n", i+1, matrixX[i]);
    }
    fflush(stdin);
    getchar();
}