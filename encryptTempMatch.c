#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

/// generator de numere pseuda-aleatoare
unsigned int generatorXor (unsigned int seed) {
    unsigned int x = seed;
    x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 5;
	return x;
}

/// structura unei imagini
typedef struct imagine {
    unsigned char *R,*G,*B,*header;
    int W, H, Wpad, size;
};

/// structura unei ferestre
typedef struct fereastra {
    char cifra, seAfiseaza;
    int iStart, iFinal, jStart, jFinal;
    double cor;
};

/// structura unui vector de ferestre
typedef struct vFerestre {
    int nr;
    struct fereastra *v;
};

/// structura unei culori
typedef struct culoare{
    unsigned char R,G,B;
};

/// salvarea si liniarizarea imaginii bmp
struct imagine salvareBitmap (char *destinatieFisier){
    struct imagine img;

    FILE *f = fopen(destinatieFisier, "rb");
    if (f == NULL){
        printf("Eroare la deschiderea fisierului!");
        return ;
    }
    else{
        /// salvam headerul
        img.header= (char *) malloc (54 * sizeof (char));
        fread(img.header, sizeof(char), 54, f);

        /// salvam dimensiunea imaginii
        fseek(f, 2, SEEK_SET);
        fread(&img.size, sizeof(unsigned int), 1, f);

        /// citim dimensiunea imaginii
        fseek(f, 18, SEEK_SET);
        fread(&img.W, sizeof(unsigned int), 1, f);
        fread(&img.H, sizeof(unsigned int), 1, f);

        /// setam dimensiunea de salvare luand in considerare pixelii de padding
        img.Wpad=(img.size-54)/img.H-3*img.W;

        /// alocam spatiu pentru pixeli
        int i, j;
        img.R=(unsigned char *) malloc (img.W*img.H*sizeof(unsigned char));
        img.G=(unsigned char *) malloc (img.W*img.H*sizeof(unsigned char));
        img.B=(unsigned char *) malloc (img.W*img.H*sizeof(unsigned char));

        /// mutam cursorul la primii pixeli si incepem citirea
        fseek(f, 54, SEEK_SET);
        for (i=0; i<img.H; ++i){
            for (j=0; j<img.W; ++j){
                fread(&(img.B[img.W*img.H-(i+1)*img.W+j]), sizeof(unsigned char), 1, f);
                fread(&(img.G[img.W*img.H-(i+1)*img.W+j]), sizeof(unsigned char), 1, f);
                fread(&(img.R[img.W*img.H-(i+1)*img.W+j]), sizeof(unsigned char), 1, f);
            }
            /// sarim peste pixelii de padding fiecare linie in parte
            fseek(f, img.Wpad, SEEK_CUR);
        }
        fclose(f);
        return img;
    }
}

/// afisarea unei imaginii bmp (salvarea ei pe harddisk)
void afisare (struct imagine img, char *destinatieSalvare){
    FILE *f;
    f = fopen(destinatieSalvare,"wb");
    if (f == NULL){
        printf("Eroare la deschiderea fisierului!");
        return ;
    }

    /// printare header
    fwrite(img.header,sizeof(char),54,f);

    /// printare pixeli
    int i, j;
    unsigned char c=0;
    for (i=0; i<img.H; ++i){
        for (j=0; j<img.W; ++j){
            fwrite(&img.B[img.W*img.H-(i+1)*img.W+j],sizeof(unsigned char),1,f);
            fwrite(&img.G[img.W*img.H-(i+1)*img.W+j],sizeof(unsigned char),1,f);
            fwrite(&img.R[img.W*img.H-(i+1)*img.W+j],sizeof(unsigned char),1,f);
        }

        /// printarea pixelilor de padding
        j=0;
        while (j!=img.Wpad){
            fwrite(&c,sizeof(unsigned char),1,f);
            ++j;
        }
    }
    fclose(f);
}

/// functie de generare a unei permutari aleatoare
unsigned int *generarePermutare (unsigned int width, unsigned int height, unsigned int *randomString){
    unsigned int i, j, aux, *perm;
    perm = (unsigned int *) malloc (width*height*sizeof(unsigned int));

    for (i=0; i<width*height; ++i){
        perm[i]=i;
    }

    for (i = width*height-1; i > 0; --i) {
        j = randomString[width*height-i] % (i + 1);
        aux = perm[j];
        perm[j] = perm[i];
        perm[i] = aux;
    }

    return perm;
}

/// functie ce calculeaza inversa unei permutari
unsigned int permutareInversa (unsigned int *permutare, unsigned int width, unsigned int height){
    unsigned int *rasp, i;
    rasp = (unsigned int *) malloc (width*height*sizeof(unsigned int));

    for (i=0; i<width*height; i++){
        rasp[permutare[i]]=i;
    }

    return rasp;
}

/// criptarea unei imagini si salvarea celei noi
void criptare (char *imagineSursa, char *imagineCriptata, char *cheie){

    FILE *sursa, *criptata, *key;
    sursa=fopen (imagineSursa,"rb");
    if (sursa == NULL){
        printf("Eroare la deschiderea fisierelor");
        return ;
    }
    key=fopen (cheie,"r");
    if (key == NULL){
        printf("Eroare la deschiderea fisierelor");
        return ;
    }

    /// citim cheile
    unsigned int R0, SV;
    fscanf(key,"%d", &R0);
    fscanf(key,"%d", &SV);

    /// citim imaginea initiala
    struct imagine imgInit;
    imgInit=salvareBitmap(imagineSursa);

    /// alocarea spatiului sirului de numere aleatoare R
    unsigned int width=imgInit.W;
    unsigned int height=imgInit.H;
    unsigned int *randomString;
    randomString = (unsigned int*) malloc (width*height*2*sizeof(unsigned int));
    randomString[0]=R0;
    int i=1;

    /// generam sirul de numere aleatoare
    while (i<(2*height*width)){
        randomString[i]=generatorXor( randomString[i-1] );
        ++i;
    }

    /// generam permutarea
    unsigned int *permutare = generarePermutare (width,height,randomString);

    /// vectori auxiliari in care vom salva permutarea pentru a o suprascrie apoi peste imaginea initiala
    unsigned char *raspunsR;
    unsigned char *raspunsG;
    unsigned char *raspunsB;
    raspunsR = (unsigned int *) malloc (width*height*sizeof(unsigned char));
    raspunsG = (unsigned int *) malloc (width*height*sizeof(unsigned char));
    raspunsB = (unsigned int *) malloc (width*height*sizeof(unsigned char));

    /// salvam in vectorii raspuns pixelii permutati
    for (i=0; i<width*height; ++i){
        raspunsB[permutare[i]]=imgInit.B[i];
        raspunsG[permutare[i]]=imgInit.G[i];
        raspunsR[permutare[i]]=imgInit.R[i];
    }

    /// copiem din vectorii raspuns in poza initiala pixelii permutati
    for (i=0; i<width*height; ++i){
        imgInit.B[i]=raspunsB[i];
        imgInit.G[i]=raspunsG[i];
        imgInit.R[i]=raspunsR[i];
    }

    /// eliberam memoria ocupata de vectorii cu pixelii permutati
    free(raspunsB);
    free(raspunsG);
    free(raspunsR);

    /// initializam SV pentru fiecare culoare si primul numar random pentru initializarea
    /// criptarii
    unsigned char SVR=0,SVG=0,SVB=0;
    for (i=0; i<=7; i++){
        SVB+=(SV>>(7-i)&1)<<(7-i);
    }
    SV=SV>>8;
    for (i=0; i<=7; i++){
        SVG+=(SV>>(7-i)&1)<<(7-i);
    }
    SV=SV>>8;
    for (i=0; i<=7; i++){
        SVR+=(SV>>(7-i)&1)<<(7-i);
    }
    unsigned char rR=0,rG=0,rB=0;
    unsigned int random = randomString[width*height];
    for (i=0; i<=7; i++){
        rB+=(random>>(7-i)&1)<<(7-i);
    }
    random=random>>8;
    for (i=0; i<=7; i++){
        rG+=(random>>(7-i)&1)<<(7-i);
    }
    random=random>>8;
    for (i=0; i<=7; i++){
        rR+=(random>>(7-i)&1)<<(7-i);
    }

    /// incepem criptarea
    imgInit.B[0]=SVB^imgInit.B[0]^rB;
    imgInit.G[0]=SVG^imgInit.G[0]^rG;
    imgInit.R[0]=SVR^imgInit.R[0]^rR;

    int j;
    for (j=1; j<width*height; ++j){
        /// calculam paramentrii random
        rR=0,rG=0,rB=0;
        random=randomString[j+width*height];
        for (i=0; i<=7; i++){
            rB+=(random>>(7-i)&1)<<(7-i);
        }
        random=random>>8;
        for (i=0; i<=7; i++){
            rG+=(random>>(7-i)&1)<<(7-i);
        }
        random=random>>8;
        for (i=0; i<=7; i++){
            rR+=(random>>(7-i)&1)<<(7-i);
        }


        imgInit.B[j]=imgInit.B[j-1]^imgInit.B[j]^rB;
        imgInit.G[j]=imgInit.G[j-1]^imgInit.G[j]^rG;
        imgInit.R[j]=imgInit.R[j-1]^imgInit.R[j]^rR;
    }

    /// eliberam vectorul cu numerele generate aleator
    free (randomString);

    /// salvam imaginea pe harddisk
    afisare(imgInit,imagineCriptata);
    fclose(sursa);
    fclose(key);
}

/// decriptarea unei imagini si salvarea celei initiale
void decriptare (char *imagineCriptata, char *imagineInitiala, char *cheie){

    FILE *sursa, *decriptata, *key;
    sursa=fopen (imagineCriptata,"rb");
    if (sursa == NULL){
        printf("Eroare la deschiderea fisierelor");
        return ;
    }
    key=fopen (cheie,"r");
    if (key == NULL){
        printf("Eroare la deschiderea fisierelor");
        return ;
    }

    /// citim cheile
    unsigned int R0, SV;
    fscanf(key,"%d", &R0);
    fscanf(key,"%d", &SV);

    /// citim imaginea criptata
    struct imagine imgCript;
    imgCript=salvareBitmap(imagineCriptata);

    /// alocarea spatiului sirului de numere aleatoare R
    unsigned int width=imgCript.W;
    unsigned int height=imgCript.H;
    unsigned int *randomString;
    randomString = (unsigned int*) malloc (width*height*2*sizeof(unsigned int));
    randomString[0]=R0;
    int i=1;

    /// generam sirul de numere aleatoare
    while (i<(2*height*width)){
        randomString[i]=generatorXor( randomString[i-1] );
        ++i;
    }

    /// generam permutarea
    unsigned int *permutare = generarePermutare (width,height,randomString);

    /// calculam inversa permutarii
    unsigned int *inversa = permutareInversa (permutare,width, height);

    /// xoram invers pixelii
    unsigned int j, random;
    unsigned char rR=0,rG=0,rB=0;
    for (j=width*height; j>0; --j){
        /// calculam paramentrii random
        rR=0,rG=0,rB=0;
        random=randomString[j+width*height];
        for (i=0; i<=7; i++){
            rB+=(random>>(7-i)&1)<<(7-i);
        }
        random=random>>8;
        for (i=0; i<=7; i++){
            rG+=(random>>(7-i)&1)<<(7-i);
        }
        random=random>>8;
        for (i=0; i<=7; i++){
            rR+=(random>>(7-i)&1)<<(7-i);
        }

        /// aplicam xorarea
        imgCript.B[j]=imgCript.B[j-1]^imgCript.B[j]^rB;
        imgCript.G[j]=imgCript.G[j-1]^imgCript.G[j]^rG;
        imgCript.R[j]=imgCript.R[j-1]^imgCript.R[j]^rR;
    }

    /// initializam SV pentru fiecare culoare si primul numar random pentru
    /// ultimul pixel de decriptat
    unsigned char SVR=0,SVG=0,SVB=0;
    for (i=0; i<=7; i++){
        SVB+=(SV>>(7-i)&1)<<(7-i);
    }
    SV=SV>>8;
    for (i=0; i<=7; i++){
        SVG+=(SV>>(7-i)&1)<<(7-i);
    }
    SV=SV>>8;
    for (i=0; i<=7; i++){
        SVR+=(SV>>(7-i)&1)<<(7-i);
    }

    random=randomString[width*height];
    rB=0, rG=0, rR=0;
    for (i=0; i<=7; i++){
        rB+=(random>>(7-i)&1)<<(7-i);
    }
    random=random>>8;
    for (i=0; i<=7; i++){
        rG+=(random>>(7-i)&1)<<(7-i);
    }
    random=random>>8;
    for (i=0; i<=7; i++){
        rR+=(random>>(7-i)&1)<<(7-i);
    }

    /// finalizam decriptarea
    imgCript.B[0]=SVB^imgCript.B[0]^rB;
    imgCript.G[0]=SVG^imgCript.G[0]^rG;
    imgCript.R[0]=SVR^imgCript.R[0]^rR;

    /// vectori auxiliari in care vom salva permutarea pentru a o suprascrie apoi peste imaginea initiala
    unsigned char *raspunsR;
    unsigned char *raspunsG;
    unsigned char *raspunsB;
    raspunsR = (unsigned int *) malloc (width*height*sizeof(unsigned char));
    raspunsG = (unsigned int *) malloc (width*height*sizeof(unsigned char));
    raspunsB = (unsigned int *) malloc (width*height*sizeof(unsigned char));

    /// salvam in vectorii raspuns pixelii permutati
    for (i=0; i<width*height; ++i){
        raspunsB[inversa[i]]=imgCript.B[i];
        raspunsG[inversa[i]]=imgCript.G[i];
        raspunsR[inversa[i]]=imgCript.R[i];
    }

    /// copiem din vectorii raspuns in poza initiala pixelii permutati
    for (i=0; i<width*height; ++i){
        imgCript.B[i]=raspunsB[i];
        imgCript.G[i]=raspunsG[i];
        imgCript.R[i]=raspunsR[i];
    }

    /// eliberam memoria ocupata de vectorii cu pixelii permutati
    free(raspunsB);
    free(raspunsG);
    free(raspunsR);

    /// eliberam vectorul cu numerele generate aleator
    free (randomString);

    /// salvam imaginea pe harddisk
    afisare(imgCript,imagineInitiala);
    fclose(sursa);
    fclose(key);
}

/// functie de calculare frecventelor canalelor de culoare
void testulChiPatrat (char *imagine){

    struct imagine img;
    img=salvareBitmap(imagine);
    printf("Chi-squared test on RGB channels for file %s: \n", imagine);

    unsigned int i;
    double fBarat, chi;
    unsigned int *frecvente;
    frecvente = (unsigned int *) calloc (256,sizeof(unsigned int));

    /// calcularea frecventelor pe canalul rosu
    chi = 0;
    for (i=0; i<img.H*img.W; ++i){
        frecvente[ img.R[i] ]++;
    }
    fBarat = img.H*img.W/256;

    for (i=0; i<256; ++i){
        chi=chi+(frecvente[i]-fBarat)*(frecvente[i]-fBarat)/fBarat;
    }
    printf ("R: %.2f\n", chi);

    /// calcularea frecventelor pe canalul verde
    memset (frecvente, 0, sizeof(unsigned int)*256);
    chi = 0;
    for (i=0; i<img.H*img.W; ++i){
        frecvente[ img.G[i] ]++;
    }
    fBarat = img.H*img.W/256;

    for (i=0; i<256; ++i){
        chi=chi+(frecvente[i]-fBarat)*(frecvente[i]-fBarat)/fBarat;
    }
    printf ("G: %.2f\n", chi);

    /// calcularea frecventelor pe canalul albastru
    memset (frecvente, 0, sizeof(unsigned int)*256);
    chi=0;
    for (i=0; i<img.H*img.W; ++i){
        frecvente[ img.B[i] ]++;
    }
    fBarat = img.H*img.W/256;

    for (i=0; i<256; ++i){
        chi=chi+(frecvente[i]-fBarat)*(frecvente[i]-fBarat)/fBarat;
    }
    printf ("B: %.2f\n", chi);

}

/// functie ce converteste o imagine la grayscale
void grayscale_image(char* nume_fisier_sursa,char* nume_fisier_destinatie) {
    FILE *fin, *fout;
    unsigned int dim_img, latime_img, inaltime_img;
    unsigned char pRGB[3], aux;

    fin = fopen(nume_fisier_sursa, "rb");
    if(fin == NULL){
        printf("nu am gasit imaginea sursa din care citesc\n");
        return;
    }

    fout = fopen(nume_fisier_destinatie, "wb+");

    fseek(fin, 2, SEEK_SET);
    fread(&dim_img, sizeof(unsigned int), 1, fin);

    fseek(fin, 18, SEEK_SET);
    fread(&latime_img, sizeof(unsigned int), 1, fin);
    fread(&inaltime_img, sizeof(unsigned int), 1, fin);
    //copiaza octet cu octet imaginea initiala in cea noua
    fseek(fin,0,SEEK_SET);
    unsigned char c;
    while(fread(&c,1,1,fin)==1){
        fwrite(&c,1,1,fout);
        fflush(fout);
    }
    fclose(fin);

    //calculam padding-ul pentru o linie
    int padding;
    if(latime_img % 4 != 0)
        padding = 4 - (3 * latime_img) % 4;
    else
        padding = 0;


    fseek(fout, 54, SEEK_SET);
    int i,j;
    for (i = 0; i < inaltime_img; i++){
        for (j = 0; j < latime_img; j++){
            //citesc culorile pixelului
            fread(pRGB, 3, 1, fout);
            //fac conversia in pixel gri
            aux = 0.299*pRGB[2] + 0.587*pRGB[1] + 0.114*pRGB[0];
            pRGB[0] = pRGB[1] = pRGB[2] = aux;
            fseek(fout, -3, SEEK_CUR);
            fwrite(pRGB, 3, 1, fout);
            fflush(fout);
        }
        fseek(fout,padding,SEEK_CUR);
    }
    fclose(fin);
    fclose(fout);
}

/// functie ce calculeaza corelatia dintre o fereastra si un sablon
double calculCorelatie (struct imagine img, struct imagine sab, int iStart, int jStart, int iFinal, int jFinal){
    int i, j;
    double suma = 0;
    double fBarat=0, sBarat=0, sigmaS=0, sigmaF=0;

    /// calculam valoarea lui sBarat
    for (i=0; i<sab.H; ++i){
        for (j=0; j<sab.W; ++j){
            sBarat+=sab.B[i*sab.W+j];
        }
    }
    sBarat=sBarat/(sab.W*sab.H);

    /// calculam valoarea lui sigmaS
    for (i=0; i<sab.H; ++i){
        for (j=0; j<sab.W; ++j){
            sigmaS=sigmaS+(sab.B[i*sab.W+j]-sBarat)*(sab.B[i*sab.W+j]-sBarat);
        }
    }
    sigmaS=sigmaS/(sab.W*sab.H-1);
    sigmaS=sqrt(sigmaS);


    /// calculam valoarea lui fBarat
    for (i=iStart; i<=iFinal; ++i){
        for (j=jStart; j<=jFinal; ++j){
            fBarat+=img.B[i*img.W+j];
        }
    }
    fBarat=fBarat/(sab.W*sab.H);

    /// calculam valoarea lui sigmaF
    for (i=iStart; i<=iFinal; ++i){
        for (j=jStart; j<=jFinal; ++j){
            sigmaF=sigmaF+(img.B[i*img.W+j]-fBarat)*(img.B[i*img.W+j]-fBarat);
        }
    }
    sigmaF=sigmaF/(sab.W*sab.H-1);
    sigmaF=sqrt(sigmaF);

    int iSab,jSab;
    for (i=iStart, iSab=0; i<=iFinal; ++i, ++iSab){
        for (j=jStart, jSab=0; j<=jFinal; ++j, ++jSab){
            suma = suma + (img.B[i*img.W+j]-fBarat)*(sab.B[iSab*sab.W+jSab]-sBarat)/(sigmaF*sigmaS);
        }
    }
    return suma/(sab.H*sab.W);
}

/// functie ce returneaza vectorul de ferestre ce au corelatia mai mare decat pragul
struct vFerestre templateMatching (struct imagine img, struct imagine sab, char cifra, double prag){
    /// variabilele auxiliare ce vor fi copiate ulterior in vectorul de ferestre returnat
    struct fereastra aux;
    struct fereastra *vAux;

    int i, j;
    double corelatie;
    struct vFerestre corespundPragului;

    int contor=0; /// numarul de ferestre ce au corelatia mai mare decat pragul
    for (i=0; i<=img.H-sab.H; ++i){
        for (j=0 ; j<=img.W-sab.W; ++j){
            corelatie = calculCorelatie (img, sab, i, j, i+sab.H-1, j+sab.W-1);
            if (corelatie > prag){ /// in cazul in care fereastra are corelatie buna incrementam contorul
                contor++;
            }
        }
    }

    vAux = (struct fereastra *) malloc (contor * sizeof (struct fereastra));/// alocam spatiul necesar salvarii ferestrelor

    contor=0;
    for (i=0; i<=img.H-sab.H; ++i){
        for (j=0 ; j<=img.W-sab.W; ++j){
            corelatie = calculCorelatie (img, sab, i, j, i+sab.H-1, j+sab.W-1);
            if (corelatie > prag){ /// in cazul in care fereastra are corelatie buna, adaugam la vectorul de ferestre
                aux.iStart=i;
                aux.iFinal=i+sab.H-1;
                aux.jStart=j;
                aux.jFinal=j+sab.W-1;
                aux.cor=corelatie;
                aux.cifra=cifra;
                aux.seAfiseaza=1;
                contor++;
                vAux[contor-1]=aux;
            }
        }
    }

    corespundPragului.nr=contor;
    corespundPragului.v=vAux;

    return corespundPragului;
}

/// functie ce primeste o imagine, o culoare si o fereastra din imagine si contureaza fereastra cu acea culoare
struct imagine conturare (struct imagine img, struct fereastra win, struct culoare col){
    int i, j;

    for (i=win.iStart; i<=win.iFinal; ++i){
        img.B[i*img.W+win.jStart]=col.B;
        img.G[i*img.W+win.jStart]=col.G;
        img.R[i*img.W+win.jStart]=col.R;
        img.B[i*img.W+win.jFinal]=col.B;
        img.G[i*img.W+win.jFinal]=col.G;
        img.R[i*img.W+win.jFinal]=col.R;
    }
    for (j=win.jStart; j<=win.jFinal; ++j){
        img.B[win.iStart*img.W+j]=col.B;
        img.G[win.iStart*img.W+j]=col.G;
        img.R[win.iStart*img.W+j]=col.R;
        img.B[win.iFinal*img.W+j]=col.B;
        img.G[win.iFinal*img.W+j]=col.G;
        img.R[win.iFinal*img.W+j]=col.R;
    }

    return img;
}

/// functia compartor al sortarii
int cmp (const void *a, const void *b){
    //struct fereastra x = *((struct fereastra *)a);
    //struct fereastra y = *((struct fereastra *)b);

    if ((*(struct fereastra *)a).cor > (*(struct fereastra *)b).cor)
        return -1;
    return 1;
}

/// functie ce returneaza minimul intre doua numere
int min (int a, int b){
    if (a<b)
        return a;
    return b;
}

/// functie ce returneaza maximul intre doua numere
int max (int a, int b){
    if (a>b)
        return a;
    return b;
}

/// functie ce calculeaza gradul de suprapunere
double suprapunere (struct fereastra a, struct fereastra b){
    int iSus,jSus,iJos,jJos;
    double rasp;

    iSus = max(a.iStart, b.iStart);
    jSus = max(a.jStart, b.jStart);
    iJos = min(a.iFinal, b.iFinal);
    jJos = min(a.jFinal, b.jFinal);

    if (iSus>iJos || jSus>jJos){
        return 0;
    }
    else{
        rasp = (iJos-iSus+1)*(jJos-jSus+1);
    }
    return rasp/( (a.iFinal-a.iStart+1)*(a.jFinal-a.jStart+1) + (b.iFinal-b.iStart+1)*(b.jFinal-b.jStart+1) - rasp);
}

/// functie ce elimina non-maximele
struct vFerestre eliminareaNonMaximelor (struct vFerestre vect){
    int i,j,contor=vect.nr;
    double x;
    struct fereastra *rasp;

    for (i=0; i<vect.nr-1; i++){
        if (vect.v[i].seAfiseaza==1){
            for (j=i+1; j<vect.nr; j++){
                if (vect.v[j].seAfiseaza==1){
                    x=suprapunere(vect.v[i], vect.v[j]);
                    if (x>0.2){
                        vect.v[j].seAfiseaza=0;
                        contor--;
                    }
                }
            }
        }
    }
    rasp = (struct ferestre *) malloc (contor * sizeof (struct fereastra));
    j=0;
    for (i=0; i<vect.nr-1; i++){
        if (vect.v[i].seAfiseaza==1){
            rasp[j]=vect.v[i];
            j++;
        }
    }
    free (vect.v);
    vect.nr=contor;
    vect.v=rasp;
    return vect;
}

int main () {
    /// partea de criptare
    char *caleImagineInitiala, *caleImagineCriptata, *caleCheieSecreta;
    caleImagineInitiala = (char *) malloc ( 250*sizeof(char) );
    caleImagineCriptata = (char *) malloc ( 250*sizeof(char) );
    caleCheieSecreta = (char *) malloc ( 250*sizeof(char) );

    printf ("Introduceti calea imaginii ce va fi criptata: ");
    scanf ("%s", caleImagineInitiala);
    printf ("Introduceti calea unde va fi salvata imaginea criptata: ");
    scanf ("%s", caleImagineCriptata);
    printf ("Introduceti calea fisierului text ce contine cheia secreta: ");
    scanf ("%s", caleCheieSecreta);

    criptare (caleImagineInitiala, caleImagineCriptata, caleCheieSecreta);

    free (caleImagineInitiala);
    free (caleImagineCriptata);
    free (caleCheieSecreta);

    /// partea de decriptare
    char *caleImagineCriptata1, *caleImagineDecriptata1, *caleCheieSecreta1;
    caleImagineCriptata1 = (char *) malloc ( 250*sizeof(char) );
    caleImagineDecriptata1 = (char *) malloc ( 250*sizeof(char) );
    caleCheieSecreta1 = (char *) malloc ( 250*sizeof(char) );

    printf ("Introduceti calea imaginii ce va fi decriptata: ");
    scanf ("%s", caleImagineCriptata1);
    printf ("Introduceti calea unde va fi salvata imaginea decriptata: ");
    scanf ("%s", caleImagineDecriptata1);
    printf ("Introduceti calea fisierului text ce contine cheia secreta: ");
    scanf ("%s", caleCheieSecreta1);

    decriptare (caleImagineCriptata1, caleImagineDecriptata1, caleCheieSecreta1);

    free (caleImagineCriptata1);
    free (caleImagineDecriptata1);
    free (caleCheieSecreta1);

    /// partea de test chi-patrat
    char *caleChi;
    caleChi = (char *) malloc ( 250*sizeof(char) );
    printf ("Introduceti calea imaginii pe care vom face testul chi-patrat: ");
    scanf ("%s", caleChi);
    testulChiPatrat(caleChi);
    free(caleChi);


    /// partea de template matching

    /// alocam spatiu pentru fiecare cale la imagine, sabloane, cat si la imaginea si sabloanele dupa aplicarea grayscale
    char *destinatie, *tempSursa, *tempSursag, *c0, *c1, *c2, *c3, *c4, *c5, *c6, *c7, *c8, *c9, *gc0, *gc1, *gc2, *gc3, *gc4, *gc5, *gc6, *gc7, *gc8, *gc9;
    destinatie = (char *) malloc ( 250*sizeof(char) );
    tempSursa = (char *) malloc ( 250*sizeof(char) );
    c0 = (char *) malloc ( 250*sizeof(char) );
    c1 = (char *) malloc ( 250*sizeof(char) );
    c2 = (char *) malloc ( 250*sizeof(char) );
    c3 = (char *) malloc ( 250*sizeof(char) );
    c4 = (char *) malloc ( 250*sizeof(char) );
    c5 = (char *) malloc ( 250*sizeof(char) );
    c6 = (char *) malloc ( 250*sizeof(char) );
    c7 = (char *) malloc ( 250*sizeof(char) );
    c8 = (char *) malloc ( 250*sizeof(char) );
    c9 = (char *) malloc ( 250*sizeof(char) );
    tempSursag = (char *) malloc ( 250*sizeof(char) );
    gc0 = (char *) malloc ( 250*sizeof(char) );
    gc1 = (char *) malloc ( 250*sizeof(char) );
    gc2 = (char *) malloc ( 250*sizeof(char) );
    gc3 = (char *) malloc ( 250*sizeof(char) );
    gc4 = (char *) malloc ( 250*sizeof(char) );
    gc5 = (char *) malloc ( 250*sizeof(char) );
    gc6 = (char *) malloc ( 250*sizeof(char) );
    gc7 = (char *) malloc ( 250*sizeof(char) );
    gc8 = (char *) malloc ( 250*sizeof(char) );
    gc9 = (char *) malloc ( 250*sizeof(char) );

    printf ("Introduceti calea unde va fi salvata imaginea dupa template matching: ");
    scanf ("%s", destinatie);
    printf ("Introduceti calea imaginii pe care vom face template matching: ");
    scanf ("%s", tempSursa);
    printf ("Introduceti calea imaginii unde se va salva versiunea grayscale a celei pe care facem template matching: ");
    scanf ("%s", tempSursag);
    printf ("Introduceti calea sabloanelor cu care vom face template matching:\n");
    scanf ("%s", c0);
    scanf ("%s", c1);
    scanf ("%s", c2);
    scanf ("%s", c3);
    scanf ("%s", c4);
    scanf ("%s", c5);
    scanf ("%s", c6);
    scanf ("%s", c7);
    scanf ("%s", c8);
    scanf ("%s", c9);
    printf ("Introduceti caile unde se vor salva versiunile grayscale ale sabloanelor cu care vom face template matching:\n");
    scanf ("%s", gc0);
    scanf ("%s", gc1);
    scanf ("%s", gc2);
    scanf ("%s", gc3);
    scanf ("%s", gc4);
    scanf ("%s", gc5);
    scanf ("%s", gc6);
    scanf ("%s", gc7);
    scanf ("%s", gc8);
    scanf ("%s", gc9);

    /// aplicam grayscale pentru toate imaginile
    grayscale_image(tempSursa, tempSursag);
    grayscale_image(c0, gc0);
    grayscale_image(c1, gc1);
    grayscale_image(c2, gc2);
    grayscale_image(c3, gc3);
    grayscale_image(c4, gc4);
    grayscale_image(c5, gc5);
    grayscale_image(c6, gc6);
    grayscale_image(c7, gc7);
    grayscale_image(c8, gc8);
    grayscale_image(c9, gc9);

    /// template matching pentru fiecare sablon in parte
    struct vFerestre rasp, detectii;
    struct imagine img, sab, imgInit;

    img=salvareBitmap(tempSursag);
    /// template matching pentru sablonul 0
    sab=salvareBitmap(gc0);
    rasp=templateMatching(img,sab,0,0.5);
    detectii.v = (struct fereastra *) malloc (rasp.nr*sizeof(struct fereastra));
    detectii.v = rasp.v;
    detectii.nr= rasp.nr;

    /// template matching pentru sablonul 1
    sab=salvareBitmap(gc1);
    rasp=templateMatching(img,sab,1,0.5);
    detectii.v = (struct fereastra *) realloc (detectii.v, (rasp.nr+detectii.nr)*sizeof(struct fereastra));
    int i;
    for (i=detectii.nr; i<detectii.nr+rasp.nr; ++i){
        detectii.v[i]=rasp.v[i-detectii.nr];
    }
    detectii.nr+=rasp.nr;

    /// template matching pentru sablonul 2
    sab=salvareBitmap(gc2);
    rasp=templateMatching(img,sab,2,0.5);
    detectii.v = (struct fereastra *) realloc (detectii.v, (rasp.nr+detectii.nr)*sizeof(struct fereastra));
    for (i=detectii.nr; i<detectii.nr+rasp.nr; ++i){
        detectii.v[i]=rasp.v[i-detectii.nr];
    }
    detectii.nr+=rasp.nr;

    /// template matching pentru sablonul 3
    sab=salvareBitmap(gc3);
    rasp=templateMatching(img,sab,3,0.5);
    detectii.v = (struct fereastra *) realloc (detectii.v, (rasp.nr+detectii.nr)*sizeof(struct fereastra));
    for (i=detectii.nr; i<detectii.nr+rasp.nr; ++i){
        detectii.v[i]=rasp.v[i-detectii.nr];
    }
    detectii.nr+=rasp.nr;

    /// template matching pentru sablonul 4
    sab=salvareBitmap(gc4);
    rasp=templateMatching(img,sab,4,0.5);
    detectii.v = (struct fereastra *) realloc (detectii.v, (rasp.nr+detectii.nr)*sizeof(struct fereastra));
    for (i=detectii.nr; i<detectii.nr+rasp.nr; ++i){
        detectii.v[i]=rasp.v[i-detectii.nr];
    }
    detectii.nr+=rasp.nr;

    /// template matching pentru sablonul 5
    sab=salvareBitmap(gc5);
    rasp=templateMatching(img,sab,5,0.5);
    detectii.v = (struct fereastra *) realloc (detectii.v, (rasp.nr+detectii.nr)*sizeof(struct fereastra));
    for (i=detectii.nr; i<detectii.nr+rasp.nr; ++i){
        detectii.v[i]=rasp.v[i-detectii.nr];
    }
    detectii.nr+=rasp.nr;

    /// template matching pentru sablonul 6
    sab=salvareBitmap(gc6);
    rasp=templateMatching(img,sab,6,0.5);
    detectii.v = (struct fereastra *) realloc (detectii.v, (rasp.nr+detectii.nr)*sizeof(struct fereastra));
    for (i=detectii.nr; i<detectii.nr+rasp.nr; ++i){
        detectii.v[i]=rasp.v[i-detectii.nr];
    }
    detectii.nr+=rasp.nr;

    /// template matching pentru sablonul 7
    sab=salvareBitmap(gc7);
    rasp=templateMatching(img,sab,7,0.5);
    detectii.v = (struct fereastra *) realloc (detectii.v, (rasp.nr+detectii.nr)*sizeof(struct fereastra));
    for (i=detectii.nr; i<detectii.nr+rasp.nr; ++i){
        detectii.v[i]=rasp.v[i-detectii.nr];
    }
    detectii.nr+=rasp.nr;

    /// template matching pentru sablonul 8
    sab=salvareBitmap(gc8);
    rasp=templateMatching(img,sab,8,0.5);
    detectii.v = (struct fereastra *) realloc (detectii.v, (rasp.nr+detectii.nr)*sizeof(struct fereastra));
    for (i=detectii.nr; i<detectii.nr+rasp.nr; ++i){
        detectii.v[i]=rasp.v[i-detectii.nr];
    }
    detectii.nr+=rasp.nr;

    /// template matching pentru sablonul 9
    sab=salvareBitmap(gc9);
    rasp=templateMatching(img,sab,9,0.5);
    detectii.v = (struct fereastra *) realloc (detectii.v, (rasp.nr+detectii.nr)*sizeof(struct fereastra));
    for (i=detectii.nr; i<detectii.nr+rasp.nr; ++i){
        detectii.v[i]=rasp.v[i-detectii.nr];
    }
    detectii.nr+=rasp.nr;
    /// in acest moment avem in vectorul detectii.v toate ferestrele care au avut prag mai mare decat 0.5

    /// etapa de sortare a vectorului de detectii in functie de gradul de corelare
    struct fereastra *deSortat=detectii.v;
    int nr = detectii.nr;
    qsort (deSortat, nr, sizeof (struct fereastra), cmp);

    /// asignarea culorilor pentru fiecare sablon
    struct culoare *vectorCulori;
    vectorCulori = (struct culoare *) malloc (10*sizeof (struct culoare));
    vectorCulori[0].R=255;
    vectorCulori[0].G=0;
    vectorCulori[0].B=0;

    vectorCulori[1].R=255;
    vectorCulori[1].G=255;
    vectorCulori[1].B=0;

    vectorCulori[2].R=0;
    vectorCulori[2].G=255;
    vectorCulori[2].B=0;

    vectorCulori[3].R=0;
    vectorCulori[3].G=255;
    vectorCulori[3].B=255;

    vectorCulori[4].R=255;
    vectorCulori[4].G=0;
    vectorCulori[4].B=255;

    vectorCulori[5].R=0;
    vectorCulori[5].G=0;
    vectorCulori[5].B=255;

    vectorCulori[6].R=192;
    vectorCulori[6].G=192;
    vectorCulori[6].B=192;

    vectorCulori[7].R=255;
    vectorCulori[7].G=140;
    vectorCulori[7].B=0;

    vectorCulori[8].R=128;
    vectorCulori[8].G=0;
    vectorCulori[8].B=128;

    vectorCulori[9].R=128;
    vectorCulori[9].G=0;
    vectorCulori[9].B=0;

    /// eliminarea nonMaximelor
    detectii=eliminareaNonMaximelor (detectii);

    /// conturarea ferestrelor ramase
    img = salvareBitmap(tempSursa);
    for (i=0; i<detectii.nr; ++i){
        conturare(img, detectii.v[i], vectorCulori[detectii.v[i].cifra]);
    }
    afisare (img, destinatie);

    free(vectorCulori);
    free(detectii.v);
    free(c0);free(c1);free(c2);free(c3);free(c4);free(c5);free(c6);free(c7);free(c8);free(c9);
    free(gc0);free(gc1);free(gc2);free(gc3);free(gc4);free(gc5);free(gc6);free(gc7);free(gc8);free(gc9);

    return 0;
}
