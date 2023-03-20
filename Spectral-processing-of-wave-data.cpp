#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <cstdlib>


#include <stdio.h>
#include <stdlib.h>


//#define TEST
//#define CopyingDataTXT



float AWGN_generator()
{/* Генерация аддитивного белого гауссовского шума с нулевым средним и стандартным отклонением, равным 1. */
 
  float temp1;
  float temp2;
  float result;
  int p;

  p = 1;

  while( p > 0 )
  {
        temp2 = ( rand() / ( (float)RAND_MAX ) ); /* функция rand() генерирует
                                                       целое число между 0 и  RAND_MAX,
                                                       которое определено в stdlib.h.
                                                   */

    if ( temp2 == 0 )
    {// temp2 >= (RAND_MAX / 2)
      p = 1;
    }// конец if
    else
    {// temp2 < (RAND_MAX / 2)
       p = -1;
    }// конец else

  }// конец while()

  temp1 = cos( ( 2.0 * (float)M_PI ) * rand() / ( (float)RAND_MAX ) );
  result = sqrt( -2.0 * log( temp2 ) ) * temp1;

  return result;        // возвращаем сгенерированный сэмпл

}


    
void FillSinusoidArray (float *const arr ,const int size) //заполнения ln
{
    int Phase=0; //фаза
    float t = 0.03; //шаг времени
    int f = 5;//частота
    float A = 2;//амплитуда сигнала
    float Interference;//белый шум

    for (int i = 0; i < size; i++)
    {   
        Interference = 0;
        // Interference = AWGN_generator();
        // if (i<300)
        // f = 1;
        // else if(i<600 && 300<i)
        // f = 3;
        // else if(i<900 && 600<i)
        // f =5;
        // else
        // f = 1;
        // arr[i]= (A*sin(2*M_PI*f*t*i+Phase)) * 100;
        arr[i]= (A*sin(2*M_PI*f*t*i+Phase));
       
    }
}

void FillnArray (float *const n ,const int size) //заполнения n
{
    for (int i = 0; i < size; i++)
    {
        n[i] = i;
    }
}

void FillnSquareArray (const float *const n, float *const nSquare ,const int size) //заполнения nSquare
{
    for (int i = 0; i < size; i++)
    {
        nSquare[i] = pow(n[i],2);
    }
}

void FilllonnArray (const float *const arr,const float *const n, float *const lonn ,const int size) //заполнения lonn
{
    for (int i = 0; i < size; i++)
    {
        lonn[i] = arr[i] * n[i];
    }
}

void FilltrnArray (const float *const n, float *const trn, const float a, const float b,const int size) //заполнения trn
{
    for (int i = 0; i < size; i++)
    {
        trn[i] = a * n[i] + b;
    }
}

void FilldlArray (float *const arr, const float *const trn, const int size) //вычисление dl
{
    for (int i = 0; i < size; i++)
    {
        arr[i] = (float)arr[i] - trn[i];
    }
}

void FillmArray (int *const m, const int size) //заполнения m
{
    for (int i = 0; i < size/2; i++)
    {
        m[i] = i;
    }
}

void ShowArray (const float *const arr, const int size)
{
   for (int i = 0; i < size; i++)
   {
       std::cout<<std::fixed<<std::setprecision(2)<<arr[i]<<"\t";
      //std::cout<<arr + i<< std::endl;
   }    
        std::cout<<std::endl;
}

void ShowArray (const int *const arr, const int size)// перегрузка функции для int
{
   for (int i = 0; i < size; i++)
   {
        std::cout<<arr[i]<<"   ";
      //std::cout<<arr + i<< std::endl;
   }    
        std::cout<<std::endl;
}

float AverageArray (const float *const arr, const int size) //среднее арифметическое массива
{
    float l_medium;
    float sum = 0;
    for (int i = 0; i < size; i++)
    {
        sum = sum + arr[i];
    }  
    l_medium = (float)sum / size;
    return l_medium;
}

float StandardDeviation( float *arr, const int size, const float l_medium) //вычисления стандартного отклонения
{
    float s = 0;
    float sum = 0;
    float divisible = 1; 
    for (int i = 0; i < size; i++)
    {
        sum = sum + pow((arr[i] - l_medium),2);
    }  
    s = sqrt((divisible / (size - 1)) * sum);
    
    return s;
}

void IteratingArray (float *const arr , const int size, const float s, const float l_medium, const int i) //для перебор массива
{
    int rightArr;
    int leftArr;
    
    for ( int j = 0; j < size; j++)
    {
        if (arr[i-j] > l_medium - 3*s && arr[i-j] < l_medium + 3*s)
        {
            leftArr = arr[i-j];
            break;            
        }
    }

    for ( int j = 0; j < size; j++)
    {
        if (arr[i+j] > l_medium - 3*s && arr[i+j] < l_medium + 3*s)
        {
            rightArr = arr[i+j];
            break;
        }
    }

    if (i != 0 && i != size-1)
    {
        arr[i] = (leftArr + rightArr) / 2;
    }
    else if (i == 0)
    {
        arr[i] = (l_medium + rightArr) / 2;
    }
    else if (i == size-1)
    {
        arr[i] = (leftArr + l_medium) / 2;
    }
}

void FilteringSeries (float *const arr , const int size, const float s, const float l_medium) //Фильтрация ряда 
{
    for (int i = 0; i < size; i++)
    {
        if (arr[i] > l_medium - 3*s && arr[i] < l_medium+ 3*s)
        {   
            arr[i] = arr[i];
        }
        else
        {
            IteratingArray(arr,size,s,l_medium,i);
        }
    }
  
}

float SignificantHeight (const float *const dln, const float dln_medium, const int size) //вычисление Hs 
{
    float Hs;
    float sum = 0;

    for (int i = 0; i < size; i++)
    {
        sum = sum + sqrt(pow(dln[i]-dln_medium,2)/(size-1));
    }

    Hs = 4.04002 * sum;

    return Hs;
}

void AutocorrFunc (float *const Km, const int *const m, const float *const dln, const int size) //Расчет автокорреляционной функции
{

    float *Ac = new float [size/2];
    int n = 0;

    for (int i = 0; i<size/2; i++)
    {
        n = 0;
        for (int z = 0; z<size-1; z++)
        {
            if (n > size)
            n = 0;
            if (z+m[i] > size)
            {
                Ac[i] = Ac[i] + dln[z] * dln[n]; 
                n++;
            }
            else
            Ac[i] = Ac[i] + dln[z] * dln[z+m[i]]; 
        }
    }

    // std::cout<<"Ac: "<< std::endl;
    // ShowArray(Ac,size/2);

    float *C = new float [size/2];

    for (int i = 0; i < (size/2); i++)
    {
        C[i] = 1 / ((float)size - m[i] -1);
    }

    // std::cout<<"C: "<< std::endl;
    // ShowArray(C,size/2);

    for (int i = 0; i < (size/2); i++)
    {
        Km[i] = Ac[i] * C[i];
    }

    delete [] C;
    delete [] Ac;
}

void SmoothingSeries (const float *const Km, float *const Kms, const int size, const int hs) //сглаживание ряда, применив механизм скользящего среднего
{
    int n = 0;
    float sum;
    for (int i=0; i<size;i++)     
    {
        sum = 0;
        for (int z=n; z<hs+n;z++)
        {
            sum += Km[z];
        }

        Kms[i] = sum / hs;
        n++;
    }

}

int FindExtr (const float *const Kms, const int size, const int hs)//выбрать все экстремумы ряда Km и рассчитать разности между минимумом и максимумом
{
    float *Extr = new float [size/2-(hs-1)];
    int *indExtr = new int [size/2-(hs-1)];//для запоминания индексов extr  массива Kms

    for (int i=0; i<size/2-(hs-1);i++) //заполняем массив 0 чтоб не было значений из кучи
    {
        Extr[i] = 0;
    }

    int Hmp_size = 0;
    int size_M;
    int maxHmp;
    int m;
    int i = 0;
    int z = 0;
    bool condition = false;

    while (z != size/2-(hs-1))//поиск экстремумов функции
    {
        if (z != 0 && z != size/2-(hs-1)-1)
        {
            if ((Kms[z] > Kms[z+1] && Kms[z] > Kms[z-1]) || (Kms[z] < Kms[z+1] && Kms[z] < Kms[z-1]))
            {
                Extr[i] = Kms[z];
                indExtr[i] = z;
                i++;
                z++;
                Hmp_size++;
            }
            else
            z++;
        }
        else if (z == 0)
        z++;
        else if (z == (size/2)-(hs-1)-1)
        z++;
    }

    // std::cout<<"================================"<<"Extr:"<<std::endl;
    // ShowArray(Extr,size/2-(hs-1));

    // std::cout<<"================================"<<"ind:"<<std::endl;
    // ShowArray(indExtr,size/2-(hs-1)); 

    if (Hmp_size == 0)//условие при котором усечение невозможно
    {
        size_M = size/2-(hs-1);
        return size_M;
    }    
    else if (Hmp_size == 1)
    {
        size_M = size/2-(hs-1);
        return size_M;
    }
    else if (Hmp_size == 2)
    {
        size_M = size/2-(hs-1);
        return size_M;
    }
    else
    {
        float *Hmp = new float [Hmp_size];

        for (int i=0; i<Hmp_size-1;i++) //заполняется массив Hmp
        {
            if (Extr[0] < 0)//Если первый экстремум в ряду  H – минимум, то значение отбрасывается.
            {
                i++;
                Hmp[i] = fabs(Extr[i+1] - Extr[i]);
            }
            else
            Hmp[i] = fabs(Extr[i+1] - Extr[i]);
        }

        // std::cout<<"================================"<<"Hmp:"<<std::endl;
        // ShowArray(Hmp,Hmp_size);   
        
        for (int i=1; i<Hmp_size-1;i++) //поиск первого от начала ряды
        {
            if (Hmp_size == 3)
            {
                if (Hmp[1] > Hmp[0])
                {
                    maxHmp = 1;
                    m = 1+1;
                    break;
                }
                else
                {
                    maxHmp = 0;
                    m = 0 + 1;
                    break;
                }
            }
            else if(Hmp[i]<Hmp[i+1] && Hmp[i]<Hmp[i-1])
            {
                maxHmp = i;
                m = i + 1 ; //так как значение m в слагаемом Extr[i+1] и есть точка усечение автокорреляционной функции 
                condition = true;
                break;
            }
            else
            condition = false;
           
        }

        // std::cout<<maxHmp<<std::endl;
        // std::cout<<m<<std::endl;
        
        if (condition == true)
        {
            for (int i=0; i<size/2-(hs-1)-1;i++)
            {
                size_M = indExtr[m];
            }
        }
        else if (condition == false)
        size_M = size/2-(hs-1);
       
        delete [] Hmp;
    }

    delete [] Extr;
    delete [] indExtr;
    return size_M;
}

void FillNewKm (const float *const Km, float *const NewKm, const int size_M,const int size_new_M )//заполнение массива NewKm + добавление 0
{
    for (int i=0; i<size_M;i++)
    {
        NewKm[i] = Km[i];
    }
    for (int i=size_M; i<size_new_M;i++)
    {
         NewKm[i] = 0;
    }

}

// void FindSk (const float *const Km, float *const Sk, const int size_M)//Расчет спектральной плотности
// {
//     float *S1 = new float [size_M];

//     float sum = 0;

//     for (int i=0; i<size_M; i++)
//     {
//         S1[i] = Km[0]/2 + ((Km[size_M-1]/2)*cos(M_PI*i));
//     }

//     std::cout<<"=========================="<<std::endl;
//     ShowArray(S1,size_M);

//     float *S2 = new float [size_M];

//     for (int i=0; i<size_M; i++)
//     {
//         sum = 0;
//         for (int z=0; z<size_M-1;z++)
//         {
//             sum = sum + (Km[z] * cos(M_PI*i*z/size_M));
//         }
//         S2[i] = sum;
//     }
//     // std::cout<<"=========================="<<std::endl;
//     // ShowArray(S2,size_M);

//     for (int i=0; i<size_M; i++)
//     {
//         Sk[i] = (S1[i] + S2[i])/M_PI;
//     }

//     delete [] S1;
//     delete [] S2;
// }


void FFTAnalysis(const float *const AVal, float *const FTvl, int const Nvl, int const Nft) //Расчет спектральной плотности
{
    int i, j, n, m, Mmax, Istp;
    float Tmpr, Tmpi, Wtmp, Theta;
    float Wpr, Wpi, Wr, Wi;
    const float TwoPi = 6.283185307179586;

    n = (int)Nvl * 2; 

    float *Tmvl = new float[n];

    for (i = 0; i < n; i+=2) 
    {
        Tmvl[i] = 0;
        Tmvl[i+1] = AVal[i/2];
    }

    i = 1; 
    j = 1;
    while (i < n) 
    {
        if (j > i)
        {
            Tmpr = Tmvl[i]; Tmvl[i] = Tmvl[j]; Tmvl[j] = Tmpr;
            Tmpr = Tmvl[i+1]; Tmvl[i+1] = Tmvl[j+1]; Tmvl[j+1] = Tmpr;
        }
        i = i + 2;
        m = Nvl;
        while ((m >= 2) && (j > m)) 
        {
            j = j - m; m = m >> 1;
        }
        j = j + m;
    }

    Mmax = 2;
    while (n > Mmax) 
    {
        Theta = -TwoPi / Mmax; Wpi = sin(Theta);
        Wtmp = sin(Theta / 2); Wpr = Wtmp * Wtmp * 2;
        Istp = Mmax * 2; Wr = 1; Wi = 0; m = 1;

        while (m < Mmax)
        {
            i = m; m = m + 2; Tmpr = Wr; Tmpi = Wi;
            Wr = Wr - Tmpr * Wpr - Tmpi * Wpi;
            Wi = Wi + Tmpr * Wpi - Tmpi * Wpr;

            while (i < n)
            {
                j = i + Mmax;
                Tmpr = Wr * Tmvl[j] - Wi * Tmvl[j-1];
                Tmpi = Wi * Tmvl[j] + Wr * Tmvl[j-1];

                Tmvl[j] = Tmvl[i] - Tmpr; Tmvl[j-1] = Tmvl[i-1] - Tmpi;
                Tmvl[i] = Tmvl[i] + Tmpr; Tmvl[i-1] = Tmvl[i-1] + Tmpi;
                i = i + Istp;
            }
        }

        Mmax = Istp;
    }

    for (i = 0; i < Nft; i++) 
    {
        j = i * 2;
        FTvl[i] = 2*sqrt(pow(Tmvl[j],2) + pow(Tmvl[j+1],2))/Nvl;
    }

    delete [] Tmvl;    
}

float FindPeakSpectrum (const float *const Sk,const int size, const float D)
{
    int indexmax;
    float value = 0;
    float Wp;

    for (int i=0; i<size;i++)
    {
        if (Sk[i] > value)
        {
            value = Sk[i];
            indexmax = i+1;
        }
    }
    Wp = M_PI*indexmax*D/size;

    // std::cout<< value << std::endl;
    // std::cout<< indexmax << std::endl;
    // std::cout<< Wp << std::endl;
    return Wp;
} 

void pop_back (float *&Array, int &size)
{
   size--;
  float *newarr = new float[size];
   for (int i = 0; i < size; i++)
      {
         newarr[i]=Array[i];
      }
   
   delete [] Array;
   Array = newarr;
}

void pop_middle (float *&Array, int &size, int const number)
{
   int l;
   float *newarr = new float[size -1];
   for (int i = 0; i < number; i++)
      {
         newarr[i]=Array[i];
      }
   for (int i = (number + 1); i < size; i++)
      {
         l = i -1;
         newarr[l]=Array[i];
      }
   delete [] Array;
   size--;
   Array = newarr;
}

void pop_front (float *&Array, int &size)
{
   size--;
   float *newarr = new float[size];
   for (int i = 0; i < size; i++)
    {
        newarr[i]=Array[i+1];
    }
   
   delete [] Array;
   Array = newarr;
}

float FinddlnsAndtime (float *const dlns, float *const time, const int size, const float t,const float D,const float Filter, float &rH_medium, float &rT_medium, float &rA_medium)
{
    float *Extr = new float [size];
    int *indExtr = new int [size];//для запоминания индексов extr  массива dlns
    float *indTime = new float [size];

    for (int i=0; i<size;i++) //заполняем массив 0 чтоб не было значений из кучи
    {
        Extr[i] = 0;
    }

     for (int i=0; i<size;i++) //заполняем массив time значениями в соответствии с шагом времени
    {
        if (i == 0)
        time[i] = 0;
        else
        time[i] = time[i-1] + t;
    }

    int size_arr_data = 0;
    int H_size;
    int T_size;
    float H_max;
    int i = 0;
    int z = 0;

    while (z != size)//поиск экстремумов функции
    {
        if (z != 0 && z != size-1)
        {
            if ((dlns[z] > dlns[z+1] && dlns[z] > dlns[z-1]) || (dlns[z] < dlns[z+1] && dlns[z] < dlns[z-1]))
            {
                Extr[i] = dlns[z];
                indTime[i] = time[z];
                indExtr[i] = z;
                i++;
                z++;
                size_arr_data++;
            }
            else
            z++;
        }
        else if (z == 0)
        z++;
        else if (z == size-1)
        z++;
    }

    // std::cout<<"================================"<<"Extr:"<<std::endl;
    // ShowArray(Extr,size);

    // std::cout<<"================================"<<"ind:"<<std::endl;
    // ShowArray(indExtr,size); 

    // std::cout<<"================================"<<"indTime:"<<std::endl;
    // ShowArray(indTime,size); 
 
    // std::cout<<"================================"<<"time:"<<std::endl;
    // ShowArray(time,size); 

     if (size_arr_data == 0 || size_arr_data == 1)//условие при котором слишком мало выборок сигнала
    {
        std::cout<<"Insufficient data"<<std::endl;
        return 0;
    }    
    if (size_arr_data == 2 && Extr[0] < 0)//условие при котором слишком мало выборок сигнала
    {
        std::cout<<"Insufficient data"<<std::endl;
        return 0;
    }
    else
    {
        H_size = size_arr_data;
        T_size = size_arr_data;
        float *H = new float [H_size];
        float *T = new float [T_size];
        float sum;

        for (int i=0; i<size_arr_data-1;i++)//заполняется массив H и T 
        {
            if (Extr[0] < 0)//Если первый экстремум в ряду  H – минимум, то значение отбрасывается.
            {
                i++;
                H[i] = fabs(Extr[i+1] - Extr[i]);
                T[i] = (indTime[i+1] - indTime[i])/D;
            }
            else
            {
                H[i] = fabs(Extr[i+1] - Extr[i]);
                T[i] = (indTime[i+1] - indTime[i])/D;
            }
        }

        // std::cout<<"================================"<<"H:"<<std::endl;
        // ShowArray(H,H_size);   
        
        // std::cout<<"================================"<<"T:"<<std::endl;
        // ShowArray(T,T_size);   
        
        for (int i=0; i<H_size+1;i++)//фильтр, отбрасывающий значения H и T, не соответствующие ветровым волнам
        {
            if (T[i] <= Filter)
            {
                if (i == 0)
                {
                    pop_front(H,H_size);
                    pop_front(T,T_size);
                    i--;
                }
                else if (i == H_size)
                {
                    pop_back(H,H_size);
                    pop_back(T,T_size);
                 
                }
                else
                {
                    pop_middle(H,H_size,i);
                    pop_middle(T,T_size,i);
                   
                }
            }
        }

        // std::cout<<"-----------------------------"<<"Hf:"<<std::endl;
        // ShowArray(H,H_size);   
        
        // std::cout<<"------------------------------"<<"Tf:"<<std::endl;
        // ShowArray(T,T_size); 

        sum = 0;
        for (int i=0; i<H_size;i++)
        {
            sum = sum + H[i];
        }
        rH_medium = sum/H_size;

        sum = 0;
        for (int i=0; i<T_size;i++)
        {
            sum = sum + T[i];
        }
        rT_medium = sum/T_size;

        rA_medium = (9.80665 * pow(rT_medium,2))/(2*M_PI);
    
        H_max = 0;
        for (int i=0; i<H_size;i++)
        {
            if (H[i]>H_max)
            H_max = H[i];
        }
    
        delete [] H;
        delete [] T;
    }
    
    delete [] Extr;
    delete [] indExtr;
    delete [] indTime;

    return H_max;
}







int main ()
{
    bool degreeof2;
    float s;
    float Hs;//значимая высота волн
    float l_medium;
    float n_medium;
    float nSquare_medium;
    float lonn_medium;
    float dln_medium;
    int size = 1024; //максимальный номер (массива)
    int size_M; //усеченная матрица
    int size_new_M;
    int check_size_M;
    float D = 1;//дискретность (Гц)
    float t = 0.03; //шаг времени с которой принимается сигнал
    float Wp;//частота пика спектра при Sk=max
    float Tp;//период пика спектра волн
    float Ap;//длина, соответствующая периоду пика спектра волн
    float H_medium;//средняя высота волн
    float T_medium;//средний период волн
    float A_medium;//средняя длина волн
    float H_max;//максимальная высота волн
    float *arr = new float [size]; 

//////////////////////////////////////////////////////////////////Фильтрация ряда +++

    FillSinusoidArray(arr,size);

    #ifdef TEST
    std::cout<<"\nln:"<<std::endl;
    ShowArray(arr,size);
    #endif


    #ifdef CopyingDataTXT
    std::ofstream signal;

    signal.open("Valies of signal",std::ios::trunc);
    signal.close();

    signal.open("Valies of signal",std::ios::app);
    if (!signal)
        {
            std::cout << "File not open" << std::endl;
            exit(1);
        }
    else
        {
            for (int i = 0; i < size; i++)
            {
                signal<<std::fixed<<std::setprecision(2)<<arr[i]<<"\n";
            }    
       
        }
    signal.close(); 
    #endif
    

    l_medium = AverageArray(arr,size);

    #ifdef TEST
    std::cout<<"\nAverageArray: "<<l_medium<<std::endl;
    #endif

    s = StandardDeviation (arr,size,l_medium);

    #ifdef TEST
    std::cout<<"\nStandardDeviation: "<<s<<std::endl;
    #endif

    FilteringSeries(arr,size,s,l_medium);

//////////////////////////////////////////////////Удаление линейного тренда +++

    #ifdef TEST
    std::cout<<"\nl`n:"<<std::endl;
    ShowArray(arr,size);
    #endif

    l_medium = AverageArray (arr,size);

    #ifdef TEST
    std::cout<<"\nl` medium: "<<l_medium<<std::endl;
    #endif

    float *n = new float [size];

    FillnArray(n,size);

    #ifdef TEST
    std::cout<<"\nn:"<<std::endl;
    ShowArray(n,size);
    #endif

    n_medium = AverageArray(n,size);

    #ifdef TEST
    std::cout<<"\nn medium: "<<n_medium<<std::endl;
    #endif    

    float *nSquare = new float [size];

    FillnSquareArray(n,nSquare,size);

    #ifdef TEST
    std::cout<<"\nn*n:"<<std::endl;
    ShowArray(nSquare,size);
    #endif

    nSquare_medium = AverageArray(nSquare,size);

    #ifdef TEST
    std::cout<<"\nn*n medium: "<<nSquare_medium<<std::endl;
    #endif    

    float *lonn = new float [size];

    FilllonnArray(arr,n,lonn,size);

    #ifdef TEST
    std::cout<<"\nl*n:"<<std::endl;
    ShowArray(lonn,size);
    #endif

    lonn_medium = AverageArray(lonn,size);

    #ifdef TEST
    std::cout<<"\nl*n medium: "<<lonn_medium<<std::endl;
    #endif

    float a = (lonn_medium - l_medium * n_medium)/(nSquare_medium - pow(n_medium,2));

    #ifdef TEST
    std::cout<<"\na: "<<a<<std::endl;
    #endif

    float b = l_medium - a * n_medium;

    #ifdef TEST
    std::cout<<"\nb: "<<b<<std::endl;
    #endif

    float *trn = new float [size];

    FilltrnArray(n,trn,a,b,size);

    #ifdef TEST
    std::cout<<"\ntrn:"<<std::endl;
    ShowArray(trn,size);
    #endif

    FilldlArray(arr,trn,size);

    #ifdef TEST
    std::cout<<"\ndln:"<<std::endl;
    ShowArray(arr,size);
    #endif

////////////////////////////////////////////////////Вычисление значимой высоты волны +++

    dln_medium = AverageArray(arr,size);

    #ifdef TEST
    std::cout<<"\ndln medium: "<<dln_medium<<std::endl;
    #endif

    Hs = SignificantHeight(arr,dln_medium,size);

    #ifdef TEST
    std::cout<<"\nHs: "<<Hs<<std::endl;
    #endif

///////////////////////////////////////////////////Расчет автокорреляционной функции +++ 
    
    int *m = new int [size/2];

    FillmArray(m,size);

    #ifdef TEST
    std::cout<<"\nm:"<<std::endl;
    ShowArray(m,size/2);
    #endif

    float *Km = new float [size/2];

    AutocorrFunc(Km,m,arr,size);

    #ifdef TEST
    std::cout<<"\nKm:"<<std::endl;
    ShowArray(Km,size/2);
    #endif

    #ifdef CopyingDataTXT
    std::ofstream list;

    list.open("Valies of Km",std::ios::trunc);
    list.close();

    list.open("Valies of Km",std::ios::app);
    if (!list)
        {
            std::cout << "File not open" << std::endl;
            exit(1);
        }
    else
        {
            for (int i = 0; i < size/2; i++)
            {
                list<<std::fixed<<std::setprecision(2)<<Km[i]<<"\n";
            }    
       
        }
    list.close(); 
    #endif

////////////////////////////////////////////Усечение автокорреляционной функции +++

    int hs = 4; //шаг осреднения 

    float *Kms = new float [size/2-(hs-1)];

    SmoothingSeries(Km,Kms,size/2-(hs-1),hs);

    #ifdef TEST
    std::cout<<"\nKms: "<<std::endl;
    ShowArray(Kms,size/2-(hs-1));
    #endif

    size_M = FindExtr(Kms,size,hs);
   
    #ifdef TEST
    std::cout<<"\nsize_M: "<<size_M<<std::endl;
    #endif

//////////////////////////////////////////Расчет спектральной плотности

    check_size_M =size_M; 
  
    if ((!(check_size_M&(check_size_M-1))) == 1) //Алгоритм проверки числа на степень 2
    degreeof2 = true;
    else
    degreeof2 = false;

    if (degreeof2 == false)//поиск нового ближайшего размера size_M кратного степени 2
    {
        size_new_M = 0;
        int i = 0;
        while (size_M > size_new_M)
        {
            i++;
            size_new_M = pow(2,i);
        }
    }
    else
    size_new_M = size_M;

    #ifdef TEST
    std::cout<<"\nsize_new_M: "<<size_new_M<<std::endl;
    #endif

    float *NewKm = new float [size_new_M];

    FillNewKm(Km,NewKm,size_M,size_new_M);

    #ifdef TEST
    std::cout<<"\nNewKm: "<<std::endl;
    ShowArray(NewKm,size_new_M);
    #endif

    float *Sk = new float [size_new_M];

    FFTAnalysis(NewKm,Sk,size_new_M,size_new_M);

    // FindSk(Km,Sk,size_M);

    #ifdef TEST
    std::cout<<"\nSk: "<<std::endl;
    ShowArray(Sk,size_M);// or size_new_M ???
    #endif

    #ifdef CopyingDataTXT
    std::ofstream list1;

    list1.open("Valies of Sk",std::ios::trunc);
    list1.close();

    list1.open("Valies of Sk",std::ios::app);
    if (!list1)
        {
            std::cout << "File not open" << std::endl;
            exit(1);
        }
    else
        {
            for (int i = 0; i < size_M; i++)
            {
                list1<<std::fixed<<std::setprecision(2)<<Sk[i]<<"\n";
            }    
       
        }
    list1.close(); 
    #endif

///////////////////////////////////////////////Расчет периода и длины волны, соответствующих пику спектра

    Wp = FindPeakSpectrum(Sk,size_M,D);
    Tp = 2*M_PI/Wp;
    Ap = (9.80665 * pow(Tp,2))/(2*M_PI);

    #ifdef TEST
    std::cout<<"Wp: "<<Wp<<std::endl;
    std::cout<<"Tp: "<<Tp<<std::endl;
    std::cout<<"Ap: "<<Ap<<std::endl;
    #endif

////////////////////////////////////////////Расчет статистических характеристик

    float Filter = 0.1; //при соответствии входного сигнала с морской волной значение должно быть 2

    #ifdef TEST
    std::cout<<"\ndln:"<<std::endl;
    ShowArray(arr,size);
    #endif

    hs = 2; //шаг осреднения 

    float *dlns = new float [size-(hs-1)];

    SmoothingSeries(arr,dlns,size-(hs-1),hs);

    #ifdef TEST
    std::cout<<"\ndlns: "<<std::endl;
    ShowArray(dlns,size-(hs-1));
    #endif

    float *time = new float [size-(hs-1)];

    float &rH_medium = H_medium;
    float &rT_medium = T_medium;
    float  &rA_medium = A_medium;

    H_max = FinddlnsAndtime(dlns,time,size-(hs-1),t,D,Filter,rH_medium,rT_medium,rA_medium);

////////////////////////////////////////////Выходные параметры

    std::cout<<std::fixed<<std::setprecision(2)<<"rH_medium: "<<H_medium<<std::endl;//средняя высота волн
    std::cout<<std::fixed<<std::setprecision(2)<<"rT_medium: "<<T_medium<<std::endl;//средний период волн
    std::cout<<std::fixed<<std::setprecision(2)<<"rA_medium: "<<A_medium<<std::endl;//средняя длина волн
    std::cout<<std::fixed<<std::setprecision(2)<<"Hs: "<<Hs<<std::endl;//значимая высота волн
    std::cout<<std::fixed<<std::setprecision(2)<<"Tp: "<<Tp<<std::endl;//период пика спектра волн
    std::cout<<std::fixed<<std::setprecision(2)<<"Ap: "<<Ap<<std::endl;//длина, соответствующая периоду пика спектра волн
    std::cout<<std::fixed<<std::setprecision(2)<<"H_max: "<<H_max<<std::endl;//максимальная высота волн

    delete [] arr;
    delete [] n;
    delete [] nSquare;
    delete [] lonn;
    delete [] trn;
    delete [] m;
    delete [] Km;
    delete [] Kms;
    delete [] NewKm;
    delete [] Sk;
    delete [] dlns;
    delete [] time;

    return 0;
}