//Anderson Zudio de Moraes
//Victor Cracel Messner
//AC em 30/03/2014

#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
using namespace std;

int primos[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};
int Fatoresprimos[100];
void Fatora(int a)
{
   for(int j=2; j<=a; j++)
   {
   int numero = j;
      while ((numero % 2) == 0)
      {
         Fatoresprimos[2]++;
         numero = numero / 2;
      }
      int i = 3;

      while (i*i <= (numero) + 1)
      {
         if ((numero % i) == 0)
         {
            Fatoresprimos[i]++;
            numero /= i;
         }
         else i = i + 2;
      }
      if (numero > 1)
         Fatoresprimos[numero]++;
   }
}
int main()
{
   int N;
   while(cin >> N && N)
   {
      memset(Fatoresprimos, 0, sizeof Fatoresprimos);
      Fatora(N);
      printf("%3ld! =",N);
      for(int j=0, c=0; j<25; j++, c++)
      {
         if(primos[j]>N)
            break;
         if(c%15==0 && c)
            printf("\n%6c",' ');
         printf("%3d",Fatoresprimos[primos[j]]);
      }
      printf("\n");
   }
   return 0;
}
