#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include "MersenneTwister.h"

using namespace std;


#define N 50
#define totalruns 100
ofstream results("results.txt");


void loadtest(vector<vector<int> >& IndX){

    ifstream datatest;
    datatest.open("testdat.txt");

    for(int i=1;i<N+1;i++){
        for(int j=1;j<N+1;j++){
            datatest>>IndX[i][j];
        }
    }

    datatest.close();
}



void seethematrix(vector<vector<int> >& IndX){

    for(int i=0;i<N+1;i++){
            //cout<<"i= "<<i<<endl;
            for(int j=0;j<N+1;j++){
               cout<<IndX[i][j]<<" ";
            }
            cout<<"\n";
        }
    cout<<"-----------------------------------"<<endl;
}



void extractit(vector<vector<int> >& IndX){

     ofstream dataout;
    dataout.open("sweepdat.txt");


    for(int i=1;i<N+1;i++){
        for(int j=1;j<N+1;j++){
           dataout<<IndX[i][j]<<"\t";
        }
        dataout<<"\n";
    }


    dataout.close();
}





void creatematrix(vector<vector<int> >& IndX,double& p){

    vector<int> zeros(N+1,0);
    fill(IndX.begin(),IndX.end(),zeros);

    MTRand random;
    int pores=p*N*N;

    int xopen=random.randInt(N-1)+1; //range [1 N], no 0
    int yopen=random.randInt(N-1)+1;

    for(int n=0;n<pores;n++){

        while(IndX[xopen][yopen]==1){
           xopen=random.randInt(N-1)+1;
           yopen=random.randInt(N-1)+1;
        }

    IndX[xopen][yopen]=1;
    }


}




void CMLT(vector<vector<int> >& IndX,vector<int>& S){


     int k=1;
    S.clear();
    S.push_back(0);//start from 1
    vector<int> L;//(N*N);
    L.push_back(0);//start from 1
    int tempL;


    for(int j=1;j<=N;j++){
        for(int i=1;i<=N;i++){

            if(IndX[i][j]!=0){

                //left
                if(IndX[i][j-1]==0){

                    //up
                    if(IndX[i-1][j]==0){

                        //L[k]=k;
                        L.push_back(k);
                        IndX[i][j]=L[k];
                        //S[L[IndX[i][j]]]++;
                        S.push_back(1);
                        k++;


                    }
                    //if up not 0
                    else{

                        IndX[i][j]=L[IndX[i-1][j]];
                        S[L[IndX[i][j]]]++;
                    }


                }
                //if left not 0
                else{

                    IndX[i][j]=L[IndX[i][j-1]];//label aristera
                    S[L[IndX[i][j]]]++;


                    if(IndX[i-1][j]!=0){//an to panw einai 0, antio

                        //L[IndX[i-1][j]]=L[IndX[i][j-1]];

                            tempL=L[IndX[i-1][j]];
                        for(int l=0;l<L.size();l++){
                            if(L[l]==tempL) {L[l]=L[IndX[i][j-1]];}
                        }

                        //merge if different only!
                        if(tempL!=L[IndX[i][j-1]]){

                        S[L[IndX[i][j-1]]]+=S[tempL];
                        S[tempL]=0;
                        }


                    }
                }


            }

        }
    }


    //second sweep
    for(int j=1;j<=N;j++){
        for(int i=1;i<=N;i++){

           if(IndX[i][j]!=0){
               IndX[i][j]=L[IndX[i][j]];
           }

        }
    }



}




bool percolation(vector<vector<int> >& IndX){

    vector<int> upper_edge = IndX[1];
    vector<int> lower_edge = IndX[N];
    vector<int> left_edge;
    for(int i=1;i<=N;i++){left_edge.push_back(IndX[i][1]);}
    vector<int> right_edge;
    for(int i=1;i<=N;i++){right_edge.push_back(IndX[i][N]);}


    for(int i=0;i<N;i++){
        if(upper_edge[i]!=0){
        for(int j=0;j<N;j++){
            if (upper_edge[i]==lower_edge[j]){return true;}
        }
        }
    }

    for(int i=0;i<N;i++){
        if(left_edge[i]!=0){
        for(int j=0;j<N;j++){
            if (left_edge[i]==right_edge[j]){return true;}
        }
        }
    }

   return false;

}




double Iav(vector<int>& S,double& p){


    int Smax=*max_element(S.begin(),S.end());
    vector<int> freq(Smax+1,0);//start from 0 include Smax
    for(int i=0;i<S.size();i++){

        freq[S[i]]++;
    }

/*
    cout<<"SMAX=  "<<Smax<<endl;
    cout<<"freq-size=  "<<freq.size()<<endl;

    for(int i=0;i<freq.size();i++){
        cout<<"i= "<<i<<"   "<<freq[i]<<endl;
    }
*/

    double sum=0;
    for(int i=0;i<freq.size();i++){
        if(freq[i]!=0){
            sum+=(double) freq[i]*i*i;
        }
    }

return sum/(p*N*N);

}




double Iavt(vector<int>& S,double& p){


    int Smax=*max_element(S.begin(),S.end());
    //int indexofSmax=distance(S.begin(),max_element(S.begin(),S.end()));
    vector<int> freq(Smax+1,0);//start from 0 include Smax
    for(int i=0;i<S.size();i++){

        freq[S[i]]++;

    }

/*
    cout<<"SMAX_new=  "<<Smax_new<<endl;
    cout<<"freq-size=  "<<freq.size()<<endl;

    for(int i=0;i<freq.size();i++){
        cout<<"i= "<<i<<"   "<<freq[i]<<endl;
    }
*/

    double sum=0;
    for(int i=0;i<freq.size();i++){
        if(freq[i]!=0){
            sum+=(double) freq[i]*i*i;
        }
    }



    sum-=(double) Smax*Smax;

return sum/(p*N*N);

}




double Pmax(vector<int>& S,double& p){

    int Smax=*max_element(S.begin(),S.end());

return Smax/(p*N*N);
}




int main()
{

    vector<int> vectors(N+1);
    vector<vector<int> > IndX(N+1,vectors);
    vector<int> S;
    int per_times;
    vector<vector<double> > Per_array(5,vector<double>());
    double sumIav=0,sumIavt=0,sumPmax=0;

    //double p=0.62;
    //loadtest(IndX);
    //seethematrix(N,IndX);


    double dp=0.1;
    for(double p=0.1;p<0.9;p=p+dp){

        if(p>0.45) dp=0.05;
        if(p>0.54) dp=0.01;
        if(p>0.63) dp=0.02;
        if(p>0.64) dp=0.05;
        if(p>0.69) dp=0.1;

        per_times=0;sumIav=0;sumIavt=0;sumPmax=0;


        for(int runs=0;runs<totalruns;runs++){


            creatematrix(IndX,p);
            //seethematrix(IndX);
            CMLT(IndX,S);



            if(percolation(IndX)==true) per_times++; //cout<<"-------------alleluia!---------------"<<endl<<endl;
            sumIav+=Iav(S,p);
            sumIavt+=Iavt(S,p);
            sumPmax+=Pmax(S,p);


        }




        Per_array[0].push_back(p);
        Per_array[1].push_back( (double) per_times/totalruns);
        Per_array[2].push_back(sumIav/totalruns);
        Per_array[3].push_back(sumIavt/totalruns);
        Per_array[4].push_back(sumPmax/totalruns);


         cout<<"p= "<<p<<" done..."<<endl;


   }//p

    cout<<endl;
    for(int q=0;q<Per_array[0].size();q++){
        cout<<"p= "<<Per_array[0][q]<<"  Percolation Mean= "<<Per_array[1][q]<<endl;
    }



    //ofstream results;
    //results.open("resultsdat300.txt");

    for(int q=0;q<Per_array[0].size();q++){
        results<<Per_array[0][q]<<"\t"<<Per_array[1][q]<<"\t"<<Per_array[2][q]<<"\t"<<Per_array[3][q]<<"\t"<<Per_array[4][q]<<endl;
    }

    //results.close();


    //seethematrix(IndX);
    //extractit(IndX);


    cout<<endl<<"results ready."<<endl;
    cout<<"If you want the rainbow, you gotta put up with the rain"<<endl;
    return 0;
}
