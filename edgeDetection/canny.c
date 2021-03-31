#include <stdio.h>
#include <math.h>
#define  PICSIZE 256
#define  MAXMASK 100

int         pic[PICSIZE][PICSIZE];
double  outpicx[PICSIZE][PICSIZE];
double  outpicy[PICSIZE][PICSIZE];
int    edgeflag[PICSIZE][PICSIZE];
double maskx[MAXMASK][MAXMASK],masky[MAXMASK][MAXMASK];
double convx[PICSIZE][PICSIZE],convy[PICSIZE][PICSIZE],
         mag[PICSIZE][PICSIZE],peaks[PICSIZE][PICSIZE],
   threshold[PICSIZE][PICSIZE];
   
int main(argc,argv)
int argc;
char **argv;
{
  int i,j,p,q,mr,centx,centy;
  double  xmask,ymask,sumx,sumy,sig,maxival,maxval,percent,HI=100,LO=30,cutoff,areaOfTops;
  FILE *fo1,*fo2,*fo3,*fp1, *fopen();
  char *foobar;
  
  argc--; argv++;
  foobar = *argv;
  fp1=fopen(foobar,"rb");
  
  argc--; argv++;
  foobar = *argv;
  fo1=fopen(foobar,"wb");
  
  fprintf(fo1, "P5\n");
  fprintf(fo1, "%d %d\n", PICSIZE, PICSIZE);
  fprintf(fo1, "255\n");
  
  argc--; argv++;
  foobar = *argv;
  sig = atof(foobar);

  argc--; argv++;
  foobar = *argv;
  fo2=fopen(foobar,"wb");

  fprintf(fo2, "P5\n");
  fprintf(fo2, "%d %d\n", PICSIZE, PICSIZE);
  fprintf(fo2, "255\n");

  argc--; argv++;
  foobar = *argv;
  fo3=fopen(foobar,"wb");

  fprintf(fo3, "P5\n");
  fprintf(fo3, "%d %d\n", PICSIZE, PICSIZE);
  fprintf(fo3, "255\n");

  argc--; argv++;
  foobar = *argv;
  percent = atof(foobar);                

  mr = (int)(sig * 3);
  centx = (MAXMASK / 2);
  centy = (MAXMASK / 2);

  for (i=0;i<256;i++)
    {
    for (j=-15;j<241;j++)
      {
      pic[i][j]  =  getc (fp1);
      }
    }
    
    // Make the masks
    for (p=-mr;p<=mr;p++)
    {
      for (q=-mr;q<=mr;q++)
      {
        xmask = q * exp(-1 * (((p * p) + (q * q)) / (2 * (sig * sig))));
        (maskx[p+centy][q+centx]) = xmask;
        ymask = p * exp(-1 * (((p * p) + (q * q)) / (2 * (sig * sig))));
        (masky[p+centy][q+centx]) = ymask;              
      }
    }
      
    // convolution
    for (i=mr;i<256-mr;i++)
    {
      for (j=mr;j<256-mr;j++)
        {
          sumx = 0;
          sumy = 0;
          for (p=-mr;p<=mr;p++)
          {
            for (q=-mr;q<=mr;q++)
            {
              sumx += pic[i+p][j+q] * maskx[p+centx][q+centy];
              sumy += pic[i+p][j+q] * masky[p+centx][q+centy];
            }
          }
          outpicx[i][j] = sumx;
          outpicy[i][j] = sumy;

          convx[i][j] = sumx;
          convy[i][j] = sumy;
        }
    }
        
    //sqrt(of squares)
    maxival = 0;
    for (i=mr;i<256-mr;i++)
    {
      for (j=mr;j<256-mr;j++)
      {
        mag[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
                                (outpicy[i][j]*outpicy[i][j])));
        if (mag[i][j] > maxival)
        {
          maxival = mag[i][j];
        }
      }
    }

    for (i=0;i<256;i++)
    {
      for (j=0;j<256;j++)
      {
        mag[i][j] = (mag[i][j] / maxival) * 255;           
        fprintf(fo1,"%c",(char)((int)(mag[i][j])));
      }
    }
    // Find the peaks: Part 2
    for(j=mr;j<256-mr;j++)
    {
      for(i=mr;i<256-mr;i++)
      {
        if((convx[i][j]) == 0.0)
        {
          convx[i][j] = .00001;
        }
        double slope = convy[i][j]/convx[i][j];
        if( (slope <= .4142)&&(slope > -.4142))
        {
          if((mag[i][j] > mag[i][j-1])&&(mag[i][j] > mag[i][j+1]))
          {
            peaks[i][j] = 255;
          }
        }
        else if( (slope <= 2.4142)&&(slope > .4142))
        {
          if((mag[i][j] > mag[i-1][j-1])&&(mag[i][j] > mag[i+1][j+1]))
          {
            peaks[i][j] = 255;
          }
        }
        else if( (slope <= -.4142)&&(slope > -2.4142))
        {
          if((mag[i][j] > mag[i+1][j-1])&&(mag[i][j] > mag[i-1][j+1]))
          {
            peaks[i][j] = 255;
          }
        }
        else
        {
          if((mag[i][j] > mag[i-1][j])&&(mag[i][j] > mag[i+1][j]))
          {
             peaks[i][j] = 255;
          }
        }
     }
  }
  for (i=0;i<256;i++)
  {
    for (j=0;j<256;j++)
    {
      fprintf(fo2,"%c",(char)((int)(peaks[i][j])));
    }
  }
  
  // Automatically get HI and LO: Part 4
  cutoff = percent*256*256;
  areaOfTops = 0;
  
  // Initialize histogram array to 0
  int histogram[256];
  for (i=0;i<256;i++)
  {
    histogram[i]=0;
  }
  for (i=0;i<256;i++)
  {
    for (j=0;j<256;j++)
    {
      histogram[(int)mag[i][j]]++;
    }
  }
  for (i=255;i>0;i--)
  {
    areaOfTops += histogram[i];
    if (areaOfTops > cutoff)
    {
      break;
    }
  }
  HI=i;
  LO=0.35*HI;
    

  // Double thresholds iterative: Part 3    
  for (int i=0;i<256;i++)
  {
    for(int j=0;j<256;j++)
    {
      if (peaks[i][j] != 0)
      {
        if (mag[i][j] > HI)
        {
          peaks[i][j] = 0;
          threshold[i][j] = 255;
        }
        else if (mag[i][j] < LO)
        {
          peaks[i][j] = threshold[i][j] = 0;
        }
      }
    }
  }
  int moretodo = 1;
  while (moretodo == 1)
  {
    moretodo = 0;
    for (int i=0;i<256;i++)
    {
      for (int j=0;j<256;j++)
      {
        if (peaks[i][j] != 0)
        {
          for (int p=-1;p<2;p++)
          {
            for(int q=-1;q<2;q++)
            {
              if (threshold[i+p][j+q] != 0)
              {
                peaks[i][j] = 0;
                threshold[i][j] = 255;
                moretodo = 1;
              }
            }
          }
        }
      }
    }
  }
  for (i=0;i<256;i++)
  {
    for (j=0;j<256;j++)
    {
      fprintf(fo3,"%c",(char)((int)(threshold[i][j])));
    }
  }    
}
