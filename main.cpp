#include<stdio.h>
#include<math.h>
#include<vector>
#include<iostream>
#include<algorithm>
#include<stdlib.h>

#define EPS 1e-6
#define RSEED 12345

using namespace std;

// ----------Handler functions------------
bool isZero(double num){
    if(abs(num) < EPS) return true;
    else return false;
}
void printM(vector<vector<double>> &W){

    int n = W.size();
    int m = W[0].size();
    for(int i = 0 ; i < n ; i++){
        for( int j = 0 ; j < m ; j++){

            printf("% 1.5f", W[i][j]);
        }
        printf("\n");
    }
    printf("\n");

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
    const double max_it = 100;
    double e = 100.0;
    vector<vector<double>> Wt;
    vector<vector<double>> Ap;
    vector<vector<double>> Ht;

    int n = A.size();
    int m = A[0].size();

    W.assign(n, vector<double>(p, 0.0)); 
    Wt.assign(p, vector<double>(n, 0.0)); 
    H.assign(p, vector<double>(m, 0.0)); 
    Ht.assign(m, vector<double>(p, 0.0)); 

    for(int i = 0 ; i < n ; i++){
        for(int j = 0 ; j < p ; j++){
            W[i][j] = double(rand()%100 + 1.0);
        }
    }

    

    for(int w = 0; w < max_it && e > s_eps ; w++){
        Ap = A;
        NormalizeCols(W, n, p);
    

        SolveSys(W, Ap, H);
       
        for(int i = 0 ; i < p ; i++){
            for(int j = 0 ; j < m ; j++) H[i][j] = max( H[i][j] , 0.0 );
        }

        

        Transpose(A, Ap);
        Transpose(H, Ht);
      

        SolveSys(Ht, Ap, Wt);        
        
        Transpose(Wt, W);

        for(int i = 0 ; i < n ; i++){
            for(int j = 0 ; j < p ; j++) W[i][j] = max( W[i][j] , 0.0 );
        }

    }

    

}

int main(){

    srand(RSEED);
    
    //Task1
    {
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
        
    }

    //Task 2
    {
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
    return 0;
}
