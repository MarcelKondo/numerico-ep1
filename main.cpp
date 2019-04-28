#include<stdio.h>
#include<math.h>
#include<vector>
#include<iostream>
#include<algorithm>
#include<stdlib.h>
#include<string>
#include<fstream>
#include<sstream>

#define EPS 1e-6
#define RSEED 12345

using namespace std;

// ----------Convenient functions------------
bool isZero(double num){
    if(abs(num) < EPS) return true;
    else return false;
}
void printM(vector<vector<double>> &W){

    int n = W.size();
    int m = W[0].size();
    for(int i = 0 ; i < n ; i++){
        for( int j = 0 ; j < m ; j++){
            if(isZero(W[i][j])) printf(" ----");
            else printf("% 1.2f", W[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

vector<vector<double>> minusM(const vector<vector<double>> &X, const vector<vector<double>> &Y){ 
    int n = X.size();
    int m = X[0].size();
    vector<vector<double>> R;
    if(n != Y.size() || m != Y[0].size()){
        printf("SUBTRACTION ERROR");
        return R;
    }
    else{
        R.assign(n, vector<double>(m, 0.0));
        for(int i = 0 ; i < n ; i++){
            for(int j = 0 ; j < m ; j++){
                R[i][j] = X[i][j] - Y[i][j];
            }
        }
        return R;
    }
}

vector<vector<double>> timesM(const vector<vector<double>> &X, const vector<vector<double>> &Y){
    int n = X.size();
    int m = Y[0].size();
    int k = X[0].size();
    vector<vector<double>> R;
    if(X[0].size() != Y.size()){
        printf("MULTIPLICATION ERROR");
        return R;
    }
    else{
        R.assign(n, vector<double>(m, 0.0));
        for(int i = 0 ; i < n ; i++){
            for(int j = 0 ; j < m ; j++){
                double aux = 0.0;
                
                for(int w = 0 ; w < k ; w++) aux += X[i][w]*Y[w][j];
                //printf("aux = %f\n", aux);
                R[i][j] = aux;
            }
        }
        //printM(R);
        return R;
    }
}


//---------------------------------------

void scCalc(double& s, double& c, double wjk, double wik){
    double tau;
    if( abs(wik) > abs(wjk) ){
        tau = (-1.0)*wjk/wik;
        c = 1/sqrt(1 + tau*tau);
        s = c*tau;
    }
    else{
        tau = (-1.0)*wik/wjk;
        s = 1/sqrt(1 + tau*tau);
        c = s*tau;
    }
    
    return;
}

void RotGivens(vector<vector<double>> &W, int n, int m, int i, int j, double s, double c){
    double aux;
    for(int r = 0 ; r < m ; r++ ){

        aux = c*W[i][r] - s*W[j][r];
        W[j][r] = s*W[i][r] + c*W[j][r];
        W[i][r] = aux;

    }
}

void SolveSys(vector<vector<double>> &W, vector<vector<double>> &b, vector<vector<double>> &x){
    int i, j, k;
    double s, c;
    int n = W.size();
    int m = W[0].size();
    
    for( k = 0 ; k < m ; k++){
        for( j = n-1 ; j >= k+1 ; j--){
            i = j-1;
            if( !isZero(W[j][k]) ){
                scCalc( s, c, W[j][k], W[i][k]);
                RotGivens( W, n, m, i, j, s, c);
                RotGivens( b, m, b[0].size(), i, j, s, c);
            }
        }
    }
    //printf("Solvesys:\n");
    //printM(W);
    //printM(b);
    x.assign(m, vector<double>(b[0].size(), 0.0));

    for(int p = 0 ; p < b[0].size() ; p++){ //solving systems in parallel
        for( k = m-1 ; k >= 0 ; k--){
            double aux = 0.0;
            
            for( j = k + 1 ;  j < m ; j++) aux += W[k][j] * x[j][p];
            
            x[k][p] = (b[k][p] - aux)/W[k][k];
           
        }
       
    }

}

void NormalizeCols(vector<vector<double>> &M, int n, int m){
    double aux;
    for(int j = 0 ; j < m ; j++){
        aux = 0.0;
        for(int i = 0 ; i < n ; i++) aux += M[i][j]*M[i][j];
        aux = sqrt(aux);
        for(int i = 0 ; i < n ; i++) M[i][j] = M[i][j]/aux;
    }
}

        
void Transpose(vector<vector<double>> &M , vector<vector<double>> &Mt){
    int n = M.size();
    int m = M[0].size();
    Mt.assign(m , vector<double>(n, 0.0));

    for(int i = 0 ; i < n ; i++){
        for(int j = 0 ; j < m ; j++){
            Mt[j][i] = M[i][j]; 
        }
    }
}

void ALS( vector<vector<double>> &A , int p, vector<vector<double>> &W, vector<vector<double>> &H){ //Alternating Least Squares method
    const double s_eps = 1e-5;
    const double max_it = 10;
    double e = -1.0; //unitialized
    vector<vector<double>> Wt;
    vector<vector<double>> Ap;
    vector<vector<double>> Ht;

    int n = A.size();
    int m = A[0].size();

    //printf("A = %d", A.size());

    W.assign(n, vector<double>(p, 0.0)); 
    Wt.assign(p, vector<double>(n, 0.0)); 
    H.assign(p, vector<double>(m, 0.0)); 
    Ht.assign(m, vector<double>(p, 0.0)); 

    for(int i = 0 ; i < n ; i++){
        for(int j = 0 ; j < p ; j++){
            W[i][j] = double(rand()%100 + 1.0);
        }
    }
    //printf("ALS 1\n");

    for(int w = 0; w < max_it && (e > s_eps || e < 0.0) ; w++){
        vector<vector<double>> last_W = W;
        Ap = A;
        NormalizeCols(W, n, p);
        //printf("ALS 3\n");
        //printf("W: %dx%d\nAp: %dx%d\nH: %dx%d\n", W.size(), W[0].size(), Ap.size(), Ap[0].size(), H.size(), H[0].size());
        SolveSys(W, Ap, H);
        //printf("ALS 2\n");
       
        for(int i = 0 ; i < p ; i++){
            for(int j = 0 ; j < m ; j++) H[i][j] = max( H[i][j] , 0.0 );
        }

        Transpose(A, Ap);
        Transpose(H, Ht);

        SolveSys(Ht, Ap, Wt);        
        
        Transpose(Wt, W);

        for(int i = 0 ; i < n ; i++){
            for(int j = 0 ; j < p ; j++){
                W[i][j] = max( W[i][j] , 0.0 );
                e = max(e, abs(W[i][j] - last_W[i][j]));
            }
        }
    }
}

void ReadMatrix(const string &file_name, vector<vector<double>> &A, int n_cols = -1, int n_lins = 784){
    ifstream file(file_name);
    if(!file){
        printf("Error opening file\n");
        return;
    }
    else{
        string s;
        A.clear();
        
        for(int i = 0 ; i < n_lins ; i++){
            int k = 0;
            int aux;
            string s;
            getline(file, s);
            stringstream ss(s);
            A.push_back(vector<double>());
            
            while(k++ != n_cols){
                //printf("k = %d ncols = %d\n", k, n_cols);
                
                ss >> aux;
                //if(aux != 0)printf("!= 0!!!!!!\n\n");
                A[i].push_back(aux/255.0);
            }
            
        }
    }
}

void TrainDigit(int digit, int n_comp, vector<vector<double>> &W,int n_cols = -1){
    vector<vector<double>> A;
    vector<vector<double>> H;
    ReadMatrix("dados_mnist/train_dig" + to_string(digit) + ".txt", A, n_cols);
    //printf("Read A: %dx%d", A.size(), A[0].size());
    //printM(A);
    //printf("foi\n");
    ALS(A, n_comp, W, H);
}

double Accuracy(){
    
}

int main(){

    srand(RSEED);
    
    //Task1
    if(1){
        int n, m;
        vector<vector<double>> W;
        vector<vector<double>> b;
        vector<vector<double>> x;

        //item a
        printf("START 1-A\n");
        n = 10;
        m = 10;
        W.assign(n, vector<double>(m, 0.0));
        b.assign(n, vector<double>(1, 1.0));
        x.assign(m, vector<double>(1, 0.0));
        for(int i = 0 ; i < n ; i++){
            for(int j = 0 ; j < m ; j++){
                if(i == j) W[i][j] = 2.0;
                else if(abs(i - j) == 1) W[i][j] = 1.0;
                else W[i][j] = 0.0;
            }
        }
        printM(W);
        SolveSys(W, b, x);
        printM(x);

        //end a
        //item b
        printf("START 1-B\n");
        n = 20;
        m = 17;
        W.assign(n, vector<double>(m, 0.0));
        b.assign(n, vector<double>(1, 1.0));
        x.assign(m, vector<double>(1, 0.0));
        for(int i = 0 ; i < n ; i++){
            for(int j = 0 ; j < m ; j++){
                if(abs(i-j) <= 4) W[i][j] = 1.0/(i+1 + j+1 -1);
                else W[i][j] = 0.0;
            }
            b[i][0] = i+1;
        }
        printM(W);
        SolveSys(W, b, x);
        printM(x);
        //end b

         //item c
        printf("START 1-C\n");
        n = 10;
        m = 10;
        W.assign(n, vector<double>(m, 0.0));
        b.assign(n, vector<double>(3, 1.0));
        x.assign(m, vector<double>(3, 0.0));
        for(int i = 0 ; i < n ; i++){
            for(int j = 0 ; j < m ; j++){
                if(i == j) W[i][j] = 2.0;
                else if(abs(i - j) == 1) W[i][j] = 1.0;
                else W[i][j] = 0.0;
            }
            b[i][1] = i+1;
            b[i][2] = 2*(i+1) - 1;
        }
        printM(W);
        SolveSys(W, b, x);
        printM(x);

        //end c
        //item d
        printf("START 1-D\n");
        n = 20;
        m = 17;
        W.assign(n, vector<double>(m, 0.0));
        b.assign(n, vector<double>(3, 1.0));
        x.assign(m, vector<double>(3, 0.0));
        for(int i = 0 ; i < n ; i++){
            for(int j = 0 ; j < m ; j++){
                if(abs(i-j) <= 4) W[i][j] = 1.0/(i+1 + j+1 -1);
                else W[i][j] = 0.0;
            }
            b[i][1] = i+1;
            b[i][2] = 2*(i+1) - 1;
        }
        printM(W);
        SolveSys(W, b, x);
        
        printM(x);
        //end d

        //item e
        printf("START 1-E\n");
        n = 5;
        m = 5;
        vector<vector<double>> Q{ vector<double>{ 2.0 , 1.0 , 1.0 , -1.0 , 1.0 },
                                  vector<double>{ 0.0 , 3.0 , 0.0 , 1.0 , 2.0 },
                                  vector<double>{ 0.0 , 0.0 , 2.0 , 2.0 , -1.0 },
                                  vector<double>{ 0.0 , 0.0 , -1.0 , 1.0 , 2.0 },
                                  vector<double>{ 0.0 , 0.0 , 0.0 , 3.0 , 1.0 }};
        b.assign(n, vector<double>(5, 1.0));
        x.assign(m, vector<double>(5, 0.0));
        SolveSys(Q, b, x);
        printM(Q);
        
    }

    //Task 2
    if(0){
        printf("START 2\n");
        vector<vector<double>> A{ vector<double>{ 3.0/10.0 , 3.0/5.0 , 0.0 },
                                  vector<double>{ 1.0/2.0  ,     0.0 , 1.0 },
                                  vector<double>{ 4.0/10.0 , 4.0/5.0 , 0.0 }};
        vector<vector<double>> W, H;
        

        ALS(A , 2, W, H);
        printf("W:\n");
        printM(W);
        printf("H:\n");
        printM(H);
                                  
    }

    //Main Task
    if(1){
        printf("START MAIN TASK:\n");
        vector<vector<double>> W[10];
        
        const int n_test = 10000; //MAX: 10000
        const int ndig_treino = 4000;
        const int p = 15;
        for(int d = 0 ; d < 10 ; d++){ //obtaining Wd for each digit
            printf("Training digit %d\n", d);
            TrainDigit(d, p, W[d], ndig_treino);
        }
        //printf("1");
        vector<vector<double>> A;
        ReadMatrix("dados_mnist/test_images.txt", A, n_test);

        //printf("2");
        vector<vector<double>> H[10];
        for(int d = 0 ; d < 10 ; d++){ //obtaining H from test_images for each digit
            vector<vector<double>> auxW = W[d];
            vector<vector<double>> auxA = A;
            printf("Obtaining H[%d]\n", d);
            SolveSys(auxW, auxA, H[d]);
        }
        //printf("3");
        int dig_pred[n_test];
        vector<double> dig_error;
        dig_error.assign(n_test, -1.0);

        printf("Predicting Digits...\n");
        for(int d = 0 ; d < 10 ; d++){  //predicting digits
            vector<vector<double>> E = minusM(A, timesM(W[d], H[d]));
            for(int k = 0 ; k < n_test ; k++){
                //printf("Predicting index %d\n", k);
                double c = 0.0;
                for(int i = 0 ; i < E.size() ; i++) c += E[i][k]*E[i][k];
                c = sqrt(c);
                if(c < dig_error[k] || dig_error[k] < 0.0){
                    dig_pred[k] = d;
                    dig_error[k] = c;
                }
            }
        }

        for(int k = 0 ; k < n_test ; k++) printf("%d ", dig_pred[k]);
        printf("\n");

        int real_labels[n_test];
        ifstream test_index("dados_mnist/test_index.txt");
        if(test_index) for(int i = 0 ; i < n_test ; i++) test_index >> real_labels[i];
        
        double accuracy = 0.0;
        double accuracy_per_digit[10];
        int digit_occur[10];
        for(int i = 0 ; i < 10 ; i++){
            accuracy_per_digit[i] = 0.0;
            digit_occur[i] = 0.0;
        }
        for(int i = 0 ; i < n_test ; i++){
            digit_occur[real_labels[i]]++;
            if(dig_pred[i] == real_labels[i]){
                accuracy += 1.0;
                accuracy_per_digit[real_labels[i]] += 1.0;
            }
        } 
        accuracy /= n_test;
        for(int i = 0 ; i < 10 ; i++){
            accuracy_per_digit[i] /= digit_occur[i];
        }

        printf("ndig_treino = %d, n_test = %d, p = %d\n", ndig_treino, n_test, p);
        printf("Accuracy for all digits: %f\n", accuracy);

        for(int i = 0 ; i < 10 ; i++){
            printf("Accuracy for digit %d: %f\n", i, accuracy_per_digit[i]);
        }
        
    }
    return 0;
}
