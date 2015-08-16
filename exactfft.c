#include "exactfft.h"

// Имя дампа (имя директории)
char* DumpName = ".NET.dump";

// Активирован режим дампа?
bool IsDumpMode = false;
//-----------------------------------------------------------------------------------------------------------------------------
double LogX(double arg, double logBase)
{
    return log(arg) / log(logBase);
}
//-----------------------------------------------------------------------------------------------------------------------------
int ToLowerPowerOf2(int arg)
{
    return (int)pow(2, (int)(LogX(arg, 2)));
}
//-----------------------------------------------------------------------------------------------------------------------------
double To_dB(double arg, double zero_db_level)
{
    return 10.0 * log(arg / zero_db_level); // log
}
//-----------------------------------------------------------------------------------------------------------------------------
double From_dB(double arg, double zero_db_level)
{
    return zero_db_level * pow(10, arg / 10.0); // exp
}
//-----------------------------------------------------------------------------------------------------------------------------
double FreqNode(double FFT_Node, double sampFreq, int N, bool isComplex)
{
    return (FFT_Node * sampFreq) / ((double)N * (isComplex ? 2.0 : 1.0));
}
//-----------------------------------------------------------------------------------------------------------------------------
double FFT_Node(double freqNode, double sampFreq, int N, bool isComplex)
{
    return (freqNode * ((double)N * (isComplex ? 2.0 : 1.0))) / sampFreq;
}
//-----------------------------------------------------------------------------------------------------------------------------
double PhaseNorm(double phase)
{
    if(phase > 0) while(phase >= M_PI)  phase -= M_2PI;
    else while(phase <= -M_PI) phase += M_2PI;

    return phase;
}
//-----------------------------------------------------------------------------------------------------------------------------
double Safe_atan2(double im, double re)
{
    return (((re < 0) ? -re : re) < FLOAT_MIN) ? 0 : atan2(im, re);
}
//-----------------------------------------------------------------------------------------------------------------------------
void fill_FFT_P(CFFT_Object* fftObj)
{
    int i, j, shift;

    // Выделяем память под вектор перестановок FFT...            
    fftObj->FFT_P = (int*)malloc(sizeof(int) * fftObj->NN);

    // Заполняем вектор изменения порядка следования данных...
    for(j = 0; j < LogX(fftObj->N, 2); ++j)
    {
        for(i = 0; i < fftObj->N; ++i)
        {
            fftObj->FFT_P[i << 1] = ((fftObj->FFT_P[i << 1] << 1) +
                    ((i >> j) & 1));
        }
    }

    shift = (fftObj->FFT_P[2] == (fftObj->N >> 1)) ? 1 : 0;
    for(i = 0; i < fftObj->NN; i += 2)
    {
        fftObj->FFT_P[i + 1] = (fftObj->FFT_P[i + 0] <<= shift) + 1;
    }
}
//-----------------------------------------------------------------------------------------------------------------------------
void fill_FFT_PP(CFFT_Object* fftObj)
{
    int i, j;

    // Выделяем память под вектор перестановок FFT...
    fftObj->FFT_PP = (int*) malloc(sizeof(int) * fftObj->NNPoly);

    // Заполняем вектор изменения порядка следования данных
    // (для полифазного FFT)...
    for(j = 0; j < LogX(fftObj->NPoly, 2); ++j)
    {
        for(i = 0; i < fftObj->NPoly; ++i)
        {
            fftObj->FFT_PP[i << 1] = ((fftObj->FFT_PP[i << 1] << 1) +
                    ((i >> j) & 1));
        }
    }

    for(i = 0; i < fftObj->NNPoly; i += 2)
    {
        fftObj->FFT_PP[i + 1] = (fftObj->FFT_PP[i + 0] <<= 1) + 1;
    }
}
//-----------------------------------------------------------------------------------------------------------------------------
void fill_FFT_TW_Cosine(CFFT_Object_t* pfftObj)
{
    // Выделяем память под взвешивающее окно...
    pfftObj->FFT_TW = (double*) malloc(fftObj->NN * sizeof(double));

    // Формирующие параметры взвешивающего косинусного окна
    double a0, a1, a2, a3, ad;

    // Характеристики окон: PS - "Peak Sidelobe" (наивысший боковой лепесток, дБ)
    switch(pfftObj->CosTW)
    {                
        case CosTW.RECTANGULAR_13dbPS:         { a0 = 1.0;       a1 = 0;         a2 = 0;         a3 = 0;         ad = 1.0;     break; }
        case CosTW.HANN_31dbPS:                { a0 = 1.0;       a1 = 1.0;       a2 = 0;         a3 = 0;         ad = 2;       break; }
        case CosTW.HAMMING_43dbPS:             { a0 = 0.54;      a1 = 0.46;      a2 = 0;         a3 = 0;         ad = 1.0;     break; }
        case CosTW.MAX_ROLLOFF_3_TERM_46dbPS:  { a0 = 0.375;     a1 = 0.5;       a2 = 0.125;     a3 = 0;         ad = 1.0;     break; }
        case CosTW.BLACKMAN_58dbPS:            { a0 = 0.42;      a1 = 0.5;       a2 = 0.08;      a3 = 0;         ad = 1.0;     break; }
        case CosTW.COMPROMISE_3_TERM_64dbPS:   { a0 = 0.40897;   a1 = 0.5;       a2 = 0.09103;   a3 = 0;         ad = 1.0;     break; }
        case CosTW.EXACT_BLACKMAN_68dbPS:      { a0 = 7938.0;    a1 = 9240.0;    a2 = 1430.0;    a3 = 0;         ad = 18608.0; break; }
        case CosTW.MIN_SIDELOBE_3_TERM_71dbPS: { a0 = 0.4243801; a1 = 0.4973406; a2 = 0.0782793; a3 = 0;         ad = 1.0;     break; }
        case CosTW.MAX_ROLLOFF_4_TERM_60dbPS:  { a0 = 10.0;      a1 = 15.0;      a2 = 6.0;       a3 = 1;         ad = 32.0;    break; }
        case CosTW.COMPROMISE1_4_TERM_82dbPS:  { a0 = 0.338946;  a1 = 0.481973;  a2 = 0.161054;  a3 = 0.018027;  ad = 1.0;     break; }
        case CosTW.COMPROMISE2_4_TERM_93dbPS:  { a0 = 0.355768;  a1 = 0.487396;  a2 = 0.144232;  a3 = 0.012604;  ad = 1.0;     break; }
        default:
        case CosTW.BLACKMAN_HARRIS_92dbPS:     { a0 = 0.35875;   a1 = 0.48829;   a2 = 0.14128;   a3 = 0.01168;   ad = 1.0;     break; }
        case CosTW.NUTTALL_93dbPS:             { a0 = 0.355768;  a1 = 0.487396;  a2 = 0.144232;  a3 = 0.012604;  ad = 1.0;     break; }
        case CosTW.BLACKMAN_NUTTALL_98dbPS:    { a0 = 0.3635819; a1 = 0.4891775; a2 = 0.1365995; a3 = 0.0106411; ad = 1.0;     break; }
        case CosTW.ROSENFIELD:                 { a0 = 0.762;     a1 = 1.0;       a2 = 0.238;     a3 = 0;         ad = a0;      break; }
    }
    double arg, wval;

    // Заполняем взвешивающее окно коэффициентами...
    for(int i = 0; i < pfftObj->N; ++i)
    {
        arg  = (2.0 * PI_M * i) / (double)pfftObj->N;
        wval = (a0 - a1 * cos(arg) + a2 * cos(2 * arg) - a3 * cos(3 * arg)) / ad;
        pfftObj->FFT_TW[(i << 1) + 1] = pfftObj->FFT_TW[(i << 1) + 0] = wval;
    }
}
//-----------------------------------------------------------------------------------------------------------------------------
double BesselI0(double arg)
{
    double numerator, denominator, z, z1, z2, z3, z4, z5,
           z6, z7, z8, z9, z10, z11, z12, z13, z_1, z_2;

    if (arg == 0.0)
    {
        return 1.0;
    }

    z = arg * arg;

    z1  = z * 0.210580722890567e-22 + 0.380715242345326e-19;
    z2  = z * z1 + 0.479440257548300e-16;
    z3  = z * z2 + 0.435125971262668e-13;
    z4  = z * z3 + 0.300931127112960e-10;
    z5  = z * z4 + 0.160224679395361e-7;
    z6  = z * z5 + 0.654858370096785e-5;
    z7  = z * z6 + 0.202591084143397e-2;
    z8  = z * z7 + 0.463076284721000e0;
    z9  = z * z8 + 0.754337328948189e2;
    z10 = z * z9 + 0.830792541809429e4;
    z11 = z * z10 + 0.571661130563785e6;
    z12 = z * z11 + 0.216415572361227e8;
    z13 = z * z12 + 0.356644482244025e9;

    numerator = z * z13 + 0.144048298227235e10;

    z_1 = z - 0.307646912682801e4;
    z_2 = z * z_1 + 0.347626332405882e7;
    denominator = z * z_2 - 0.144048298227235e10;

    return -numerator / denominator;
}
//-----------------------------------------------------------------------------------------------------------------------------
void fill_FFT_TW_Kaiser(CFFT_Object_t* pfftObj)
{
    int i, j;
    double norm, arg, w;

    // Выделяем память под вектор перестановок FFT...
    pfftObj->FFT_TW = (double*) malloc(pfftObj->NN * sizeof(double));

    // Нормирующий коэффициент окна Кайзера
    norm = BesselI0(pfftObj->Beta);

    // Заполняем взвешивающее окно...
    for(i = 1; i <= (pfftObj->N >> 1); ++i)
    {
        // arg = Beta * sqrt(1-(((2*(i-1))/(N-1))-1)^2);
        arg = pfftObj->Beta *
            sqrt(
                    1 - pow(
                        (
                         (double)((i - 1) << 1)
                         /
                         (double)(pfftObj->N - 1)
                        ) - 1
                        , 2)
                );

        w = BesselI0(arg) / norm;

        j = i - 1; // Приводим индекс от базы "1" к базе "0"
        pfftObj->FFT_TW[(j << 1) + 0] = w; // left re
        pfftObj->FFT_TW[(j << 1) + 1] = w; // left im
        pfftObj->FFT_TW[(fftObj->NN - 2) - (j << 1) + 0] = w; // right re
        pfftObj->FFT_TW[(fftObj->NN - 2) - (j << 1) + 1] = w; // right im
    }
}
//-----------------------------------------------------------------------------------------------------------------------------
bool CFFT_Inspector(CFFT_Object_t* pfftObj)
{
    // Размер окна не может быть меньше предельно-допустимого!
    if((pfftObj->N     < MIN_FRAME_WIDTH) ||
       (pfftObj->NPoly < MIN_FRAME_WIDTH) ||
       (pfftObj->Beta  > MAX_KAISER_BETA))
    {
        return false;
    }

    return true;
}
//-----------------------------------------------------------------------------------------------------------------------------
void CFFT_Destructor(CFFT_Object* pfftObj)
{
    pfftObj->FFT_P  = NULL;
    pfftObj->FFT_PP = NULL;
    pfftObj->FFT_TW = NULL;
}
//-----------------------------------------------------------------------------------------------------------------------------
CFFT_Object_t* CFFT_Init(int frameWidth, CosTW cosTW, double beta, int polyDiv2,
        int windowStep, bool isComplex)
{
    // Объект-результат
    CFFT_Object* pfftObj = malloc(sizeof(CFFT_Object_t));

    // Заполнение полей объекта
    pfftObj->N = ToLowerPowerOf2(frameWidth); // Размер кадра FFT
    pfftObj->NN = pfftObj->N << 1;              // Кол-во точек (re + im)
    pfftObj->NPoly = pfftObj->N >> polyDiv2;    // Размер полифазного кадра FFT
    pfftObj->NNPoly = pfftObj->NPoly << 1;      // Кол-во точек полифазного FFT
    pfftObj->CosTW = cosTW;                   // Тип косинусного взвешивающего окна
    pfftObj->Beta = beta;                     // Форм-ий коэфф. окна Кайзера
    pfftObj->PolyDiv = 1 << polyDiv2;         // Полифазный делитель
    pfftObj->WindowStep = windowStep;         // Шаг окна FFT
    pfftObj->IsComplex = isComplex;           // Используется комплексный вход? (бывает COMPLEX и L+R)

    fill_FFT_P(pfftObj);  // Вектор изменения порядка след. данных перед FFT
    fill_FFT_PP(pfftObj); // Вектор изменения порядка... (для полифазного FFT)

    if(pfftObj->CosTW == CosTW.NONE) //...если не задано взвешивающее окно косинусного типа
    {
        fill_FFT_TW_Kaiser(pfftObj); // Взвешивающее окно Кайзера

    } else
    {
        fill_FFT_TW_Cosine(pfftObj); // Косинусное взвешивающее окно
    }

    // Обрабатываем ситуацию со сбросом дампа...
    /* if(IsDumpMode)
    {
        Directory.CreateDirectory(DumpName);

        DebugHelper.WriteInts(DumpName,    "FFT_P.int32",   pfftObj.FFT_P);
        DebugHelper.WriteInts(DumpName,    "FFT_PP.int32",  pfftObj.FFT_PP);
        DebugHelper.WriteDoubles(DumpName, "FFT_TW.double", pfftObj.FFT_TW);
    } */

    // Если некоторые параметры не соответствуют норме...
    if(!CFFT_Inspector(pfftObj))
    {
        //...- убираем объект.
        CFFT_Destructor(pfftObj);
    }

    // Возвращаем объект "FFT"
    return pfftObj;
}
//-----------------------------------------------------------------------------------------------------------------------------
CFFT_Object* CFFT_Constructor_Cosine(int frameWidth, CosTW cosTW, int polyDiv2, int windowStep, bool isComplex)
{
    // Возвращаем объект "FFT"
    return CFFT_Init(frameWidth, cosTW, MAX_KAISER_BETA, polyDiv2, windowStep, isComplex);
}
//-----------------------------------------------------------------------------------------------------------------------------
CFFT_Object* CFFT_Constructor_Kaiser(int frameWidth, double beta, int polyDiv2, int windowStep, bool isComplex)
{
    // Возвращаем объект "FFT"
    return CFFT_Init(frameWidth, CosTW.NONE, beta, polyDiv2, windowStep, isComplex);
}
//-----------------------------------------------------------------------------------------------------------------------------
void CFFT_Process(double* FFT_S, int FFT_S_Offset, double* FFT_T, bool useTaperWindow, bool recoverAfterTaperWindow,
        bool useNorm, bool direction, bool usePolyphase, CFFT_Object* pfftObj)
{
    int i, j, mmax, isteps, NN, istep, ii, m, jj;
    double isign, theta, sin05Th, wpr, wpi, wr, wi, tempr, tempi, wtemp;

    // Использование взвешивающего окна допустимо
    // только при прямом преобразовании
    if(direction && useTaperWindow)
    {
        if(!usePolyphase)
        {
            // Обычное FFT
            // только тогда, когда оно активно
            for(i = 0; i < pfftObj->NN; ++i)
            {
                FFT_T[i] = pfftObj->FFT_TW[pfftObj->FFT_P[i]] *
                    FFT_S[pfftObj->FFT_P[i] + FFT_S_Offset];
            }
        }
        else
        {
            // Полифазное FFT
            // предполагает применение только на прямом проходе
            for(i = 0; i < pfftObj->NNPoly; ++i)
            {
                FFT_T[i] = 0;

                // Накапливаем сумму текущей точки (в соответствии
                // с количеством сегментов)
                for(j = 0; j < pfftObj->PolyDiv; ++j)
                {
                    FFT_T[i] += pfftObj->FFT_TW[pfftObj->FFT_PP[i] +
                        (j * pfftObj->NNPoly)] * FFT_S[pfftObj->FFT_PP[i] +
                        (j * pfftObj->NNPoly)  + FFT_S_Offset];
                }
            }
        }
    }
    else
    {
        // Обратный проход или прямой...
        // но без взвешивающего окна
        for(i = 0; i < pfftObj->NN; ++i)
        {
            FFT_T[i] = FFT_S[pfftObj->FFT_P[i] + FFT_S_Offset];
        }
    }

    // Нормализация коэффициентов производится при её выборе и только на
    // прямом проходе алгоритма (или если ситуация 100% симметрична)
    if((!direction) && (!useNorm))
    {
        for(i = 0; i < pfftObj->NNPoly; ++i)
        {
            FFT_T[i] /= pfftObj->N;
        }
    }

    // FFT Routine
    isign  = direction ? -1 : 1;
    mmax   = 2;
    isteps = 1;
    NN = usePolyphase ? pfftObj->NNPoly : pfftObj->NN;
    while(NN > mmax)
    {
        isteps++;
        istep   = mmax << 1;
        theta   = isign * ((2 * M_PI) / mmax);
        sin05Th = sin(0.5 * theta);
        wpr = -(2.0 * (sin05Th * sin05Th));
        wpi = sin(theta);
        wr  = 1.0;
        wi  = 0.0;

        for(ii = 1; ii <= (mmax >> 1); ++ii)
        {
            m = (ii << 1) - 1;
            for(jj = 0; jj <= ((NN - m) >> isteps); ++jj)
            {
                i = m + (jj << isteps);
                j = i + mmax;
                tempr = wr * FFT_T[j - 1] - wi * FFT_T[j];
                tempi = wi * FFT_T[j - 1] + wr * FFT_T[j];
                FFT_T[j - 1]  = FFT_T[i - 1] - tempr;
                FFT_T[j - 0]  = FFT_T[i - 0] - tempi;
                FFT_T[i - 1] += tempr;
                FFT_T[i - 0] += tempi;
            }
            wtemp = wr;
            wr = wr * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }

    // Нормализация коэффициентов производится при её выборе и только
    // на прямом проходе алгоритма (или если ситуация 100% симметрична)
    if(direction && useNorm)
    {
        for(i = 0; i < NN; ++i)
        {
            FFT_T[i] /= pfftObj->N;
        }
    }

    // Аннигилируем взвешивающее окно (если оно равно нулю в некоторых точках
    // - результат восстановления неизвестен и полагаем его равным нулю)
    if((!direction) && useTaperWindow && recoverAfterTaperWindow)
    {
        for(i = 0; i < pfftObj->NN; ++i)
        {
            FFT_T[i] = ((pfftObj->FFT_TW[i] == 0) ?
                    0 : (FFT_T[i] / pfftObj->FFT_TW[i]));
        }
    }
}
//-----------------------------------------------------------------------------------------------------------------------------
void CFFT_Explore(double* FFT_T, double* MagL, double* MagR, double* ACH, double* ArgL, double* ArgR,
        double* PhaseLR, bool usePolyphase, CFFT_Object* pfftObj)
{
    int i, N;
    double magL, magR, FFT_T_i_Re, FFT_T_i_Im, FFT_T_N_i_Re, FFT_T_N_i_Im,
           lx, ly, rx, ry, argL, argR;

    // Определяем размерность FFT
    N = usePolyphase ? pfftObj->NPoly : pfftObj->N;

    // Вычисляем значения искомых величин для начальной точки
    // ("нулевая" гармоника)
    magL = FFT_T[0];
    magR = FFT_T[1];
    if(MagL    != null) MagL[0]    = magL;
    if(MagR    != null) MagR[0]    = magR;
    if(ACH     != null) ACH[0]     = magR / ((magL == 0) ? FLOAT_MIN : magL);
    if(ArgL    != null) ArgL[0]    = M_PI;
    if(ArgR    != null) ArgR[0]    = M_PI;
    if(PhaseLR != null) PhaseLR[0] = 0;

    // Работа с гармоническими точками
    for(i = 1; i < (N >> 1); ++i)
    {
        FFT_T_i_Re   = FFT_T[(i << 1) + 0];
        FFT_T_i_Im   = FFT_T[(i << 1) + 1];
        FFT_T_N_i_Re = FFT_T[((N - i) << 1) + 0];
        FFT_T_N_i_Im = FFT_T[((N - i) << 1) + 1];

        lx = FFT_T_i_Re   + FFT_T_N_i_Re;
        ly = FFT_T_i_Im   - FFT_T_N_i_Im;
        rx = FFT_T_i_Im   + FFT_T_N_i_Im;
        ry = FFT_T_N_i_Re - FFT_T_i_Re;

        magL = sqrt((lx * lx) + (ly * ly)) * 0.5;
        magR = sqrt((rx * rx) + (ry * ry)) * 0.5;
        argL = Safe_atan2(ly, lx);
        argR = Safe_atan2(ry, rx);

        if(MagL    != null) MagL[i]    = magL;
        if(MagR    != null) MagR[i]    = magR;
        if(ACH     != null) ACH[i]     = magR / ((magL == 0) ? FLOAT_MIN : magL);
        if(ArgL    != null) ArgL[i]    = argL;
        if(ArgR    != null) ArgR[i]    = argR;
        if(PhaseLR != null) PhaseLR[i] = PhaseNorm(argR - argL);
    }
}
//-----------------------------------------------------------------------------------------------------------------------------
void CFFT_ComplexExplore(double* FFT_T, double* Mag, double* Arg,
        bool usePolyphase, bool isMirror, CFFT_Object* pfftObj)
{
    int i, N, N_2;
    double magL, magR, FFT_T_i_Re, FFT_T_i_Im, FFT_T_N_i_Re, FFT_T_N_i_Im,
           lx, ly, rx, ry, argL, argR;

    // Определяем размерность FFT
    N   = usePolyphase ? pfftObj->NPoly : pfftObj->N;
    N_2 = N >> 1;

    // Вычисляем значения искомых величин для начальной точки
    // ("нулевая" гармоника)
    lx = FFT_T[0];
    ly = FFT_T[1];
    rx = FFT_T[N + 0];
    ry = FFT_T[N + 1];

    if(Mag != null)
    {
        magL = sqrt((lx * lx) + (ly * ly));
        magR = sqrt((rx * rx) + (ry * ry));

        if(!isMirror)
        {
            Mag[0]   = magL;
            Mag[N_2] = magR;

        } else
        {
            Mag[N_2 - 1] = magL;
            Mag[N   - 1] = magR;
        }
    }

    if(Arg != null)
    {
        argL = Safe_atan2(ly, lx);
        argR = Safe_atan2(ry, rx);

        if(!isMirror)
        {
            Arg[0]   = argL;
            Arg[N_2] = argR;

        } else
        {
            Arg[N_2 - 1] = argL;
            Arg[N   - 1] = argR;
        }
    }

    // Работа с гармоническими точками
    for(i = 1; i < N_2; ++i)
    {
        FFT_T_i_Re   = FFT_T[(i << 1) + 0];
        FFT_T_i_Im   = FFT_T[(i << 1) + 1];
        FFT_T_N_i_Re = FFT_T[((N - i) << 1) + 0];
        FFT_T_N_i_Im = FFT_T[((N - i) << 1) + 1];

        lx = FFT_T_i_Re;
        ly = FFT_T_i_Im;
        rx = FFT_T_N_i_Re;
        ry = FFT_T_N_i_Im;

        if(Mag != null)
        {
            magL = sqrt((lx * lx) + (ly * ly));
            magR = sqrt((rx * rx) + (ry * ry));

            if(!isMirror)
            {
                Mag[i]     = magL;
                Mag[N - i] = magR;

            } else
            {
                Mag[N_2 - i - 1] = magL;
                Mag[N_2 + i - 1] = magR;
            }
        }

        if(Arg != null)
        {
            argL = Safe_atan2(ly, lx);
            argR = Safe_atan2(ry, rx);

            if(!isMirror)
            {
                Arg[i]     = argL;
                Arg[N - i] = argR;

            } else
            {
                Arg[N_2 - i - 1] = argL;
                Arg[N_2 + i - 1] = argR;
            }
        }                
    }
}
//-----------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Перевод значений массива double в форму dB
/// </summary>
/// <param name="data"> Данные для обработки. </param>
/// <param name="zero_db_level"> Уровень "нулевого" уровня. </param>
/// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
void dB_Scale(double* Mag, double zero_db_level, CFFT_Object* pfftObj)
{
    int i;
    for(i = 0; i < (fftObj.N * ((pfftObj->IsComplex ? 2 : 1))) >> 1; ++i)
    {
        Mag[i] = 10.0 * log(Mag[i] / zero_db_level); // log
    }
}
//-----------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Самотестирование внутренней точности в цикле прямого-обратного
/// преобразования на пользовательских данных
/// </summary>
/// <param name="FFT_S"> Вектор входных данных ("левый" и "правый" каналы
/// - чет./нечет.) </param>
/// <param name="ACH_Difference"> Коэффициент превосходства правого канала
/// по уровню в ходе проводимого теста. </param>
/// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
/// <returns> Структура "Результат тестирования внутренней точности
/// прямого-обратного комплексного преобразования Фурье". </returns>
CFFT_SelfTestResult SelfTest_S(double* FFT_S, double ACH_Difference, CFFT_Object* pfftObj)
{
    double* FFT_S_backward, FFT_T, MagL, MagR, ACH, ArgL, ArgR, PhaseLR;
    int N2, FFT_S_Offset, i;
    bool useTaperWindow, recoverAfterTaperWindow, useNorm, direction, usePolyphase;
    double maxDiff, FFT_T_i_Re, FFT_T_i_Im, FFT_T_N_i_Re, FFT_T_N_i_Im,
           lx, ly, rx, ry, lx_, ly_, rx_, ry_, currentDiff;
    long startCounter, CFFT_Process_counter, CFFT_Explore_counter, timerFrequency;
    int N_iters = 10000;

    // Cтруктура "Результат тестирования внутренней точности
    // прямого-обратного комплексного преобразования Фурье"
    CFFT_SelfTestResult_t* selfTestResult = (CFFT_SelfTestResult_t*) malloc(sizeof(CFFT_SelfTestResult_t));

    // Массив исходных данных - для заполнения на обратном ходе FFT
    FFT_S_backward = (double*) malloc(pfftObj->NN * sizeof(double));

    // Целевой массив
    FFT_T = (double*) malloc(pfftObj->NN * sizeof(double));

    // (Количество точек FFT / 2) - количество гармоник вместе с нулевой
    N2 = pfftObj->N >> 1;

    // Массивы результатов Фурье-анализа
    MagL    = (double*) malloc(N2 * sizeof(double));
    MagR    = (double*) malloc(N2 * sizeof(double));
    ACH     = (double*) malloc(N2 * sizeof(double));
    ArgL    = (double*) malloc(N2 * sizeof(double));
    ArgR    = (double*) malloc(N2 * sizeof(double));
    PhaseLR = (double*) malloc(N2 * sizeof(double));

    // Не используем взвешивающее окно, но работаем
    // с нормализацией - направление прямое
    useTaperWindow = FALSE;
    FFT_S_Offset   = 0;
    recoverAfterTaperWindow = FALSE;
    useNorm      = TRUE;
    direction    = TRUE;
    usePolyphase = FALSE;
    CFFT_Process(FFT_S, FFT_S_Offset, FFT_T, useTaperWindow, recoverAfterTaperWindow, useNorm, direction,
            usePolyphase, pfftObj);

    // Пользуясь результатами прямого преобразования извлекаем
    // все возможные величины
    CFFT_Explore(FFT_T, MagL, MagR, ACH, ArgL, ArgR, PhaseLR,
            usePolyphase, pfftObj);

    // Проверяем правильность расчета магнитуд и фаз - пытаемся получить
    // исходные комплексные числа
    maxDiff = 0;
    for(i = 1; i < N2; ++i)
    {
        FFT_T_i_Re   = FFT_T[(i << 1) + 0];
        FFT_T_i_Im   = FFT_T[(i << 1) + 1];
        FFT_T_N_i_Re = FFT_T[((pfftObj->N - i) << 1) + 0];
        FFT_T_N_i_Im = FFT_T[((pfftObj->N - i) << 1) + 1];

        lx = FFT_T_i_Re   + FFT_T_N_i_Re;
        ly = FFT_T_i_Im   - FFT_T_N_i_Im;
        rx = FFT_T_i_Im   + FFT_T_N_i_Im;
        ry = FFT_T_N_i_Re - FFT_T_i_Re;

        lx_ = 2 * MagL[i] * cos(ArgL[i]);
        ly_ = 2 * MagL[i] * sin(ArgL[i]);
        rx_ = 2 * MagR[i] * cos(ArgR[i]);
        ry_ = 2 * MagR[i] * sin(ArgR[i]);

        currentDiff = abs(lx - lx_);
        maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;

        currentDiff = abs(ly - ly_);
        maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;

        currentDiff = abs(rx - rx_);
        maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;

        currentDiff = abs(ry - ry_);
        maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;
    }

    // Сохраняем максимальную невязку перехода из алгебраической формы
    // в показательную и обратно
    selfTestResult->MaxDiff_ALG_to_EXP_to_ALG = maxDiff;

    // Не используем взвешивающее окно, но работаем
    // с нормализацией - направление обратное
    useTaperWindow = false;
    FFT_S_Offset   = 0;
    recoverAfterTaperWindow = false;
    useNorm      = true;
    direction    = false;
    usePolyphase = false;
    CFFT_Process(FFT_T, FFT_S_Offset, FFT_S_backward, useTaperWindow,
            recoverAfterTaperWindow, useNorm, direction,
            usePolyphase, pfftObj);

    maxDiff = 0;
    for(i = 0; i < pfftObj->N; ++i)
    {
        currentDiff = abs(FFT_S_backward[i] - FFT_S[i]);
        maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;
    }

    // Сохраняем максимальную невязку после прямого-обратного преобразования
    // (без взвешивающего окна)
    selfTestResult->MaxDiff_FORWARD_BACKWARD = maxDiff;

    // Используем взвешивающее окно, работаем
    // с нормализацией - направление прямое
    useTaperWindow = true;
    FFT_S_Offset   = 0;
    recoverAfterTaperWindow = false;
    useNorm      = true;
    direction    = true;
    usePolyphase = false;
    CFFT_Process(FFT_S, FFT_S_Offset, FFT_T, useTaperWindow,
            recoverAfterTaperWindow, useNorm, direction,
            usePolyphase, pfftObj);

    // Используем аннигиляцию взвешивающего окна, работаем
    // с нормализацией - направление обратное
    useTaperWindow = true;
    FFT_S_Offset   = 0;
    recoverAfterTaperWindow = true;
    useNorm      = true;
    direction    = false;
    usePolyphase = false;
    CFFT_Process(FFT_T, FFT_S_Offset, FFT_S_backward, useTaperWindow,
            recoverAfterTaperWindow, useNorm, direction,
            usePolyphase, pfftObj);

    maxDiff = 0;
    for(i = (pfftObj->NN / 2); i <= ((pfftObj->NN * 3) / 4); ++i)
    {
        currentDiff = abs(FFT_S_backward[i] - FFT_S[i]);
        maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;
    }

    // Сохраняем максимальную невязку после прямого-обратного
    // преобразования (c аннигиляцией взвешивающего окна)
    selfTestResult->MaxDiff_FORWARD_BACKWARD_AntiTW = maxDiff;

    maxDiff = 0;
    for(i = 0; i < N2; ++i)
    {
        currentDiff = abs(ACH[i] - ACH_Difference);
        maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;
    }

    // Сохраняем максимальную невязку по расчету заданной АЧХ
    selfTestResult->MaxDiff_ACH = maxDiff;

    maxDiff = 0;
    for(i = 0; i < N2; ++i)
    {
        currentDiff = abs(PhaseLR[i]);
        maxDiff     = (maxDiff < currentDiff) ? currentDiff : maxDiff;
    }

    // Сохраняем максимальную невязку по расчету разности хода фаз каналов
    selfTestResult->MaxDiff_PhaseLR = maxDiff;

    // Performance Test
    timerFrequency = 10000000;
    CFFT_Process_counter = CFFT_Explore_counter = 0;

    // Не используем взвешивающее окно, но работаем
    // с нормализацией - направление прямое
    useTaperWindow = false;
    FFT_S_Offset   = 0;
    recoverAfterTaperWindow = false;
    useNorm      = false;
    direction    = true;
    usePolyphase = false;

    // CFFT_Process_time
    startCounter = 0;
    startCounter = DateTime.Now.Ticks;
    for(i = 0; i < N_iters; ++i)
    {
        CFFT_Process(FFT_S, FFT_S_Offset, FFT_T, useTaperWindow,
                recoverAfterTaperWindow, useNorm, direction,
                usePolyphase, pfftObj);
    }
    CFFT_Process_counter  = DateTime.Now.Ticks;
    CFFT_Process_counter -= startCounter;
    selfTestResult->CFFT_Process_time  = (double)CFFT_Process_counter / (double)timerFrequency;
    selfTestResult->CFFT_Process_time /= (double)N_iters;

    // CFFT_Explore_time
    startCounter = 0;
    startCounter = DateTime.Now.Ticks;
    for(i = 0; i < N_iters; ++i)
    {
        CFFT_Explore(FFT_T, MagL, MagR, ACH, ArgL, ArgR, PhaseLR,
                usePolyphase, pfftObj);
    }
    CFFT_Explore_counter  = DateTime.Now.Ticks;
    CFFT_Explore_counter -= startCounter;
    selfTestResult->CFFT_Explore_time  = (double)CFFT_Explore_counter / (double)timerFrequency;
    selfTestResult->CFFT_Explore_time /= (double)N_iters;

    // Высвобождаем ресурсы динамической памяти
    FFT_S = null;
    FFT_S_backward = null;
    FFT_T   = null;
    MagL    = null;
    MagR    = null;
    ACH     = null;
    ArgL    = null;
    ArgR    = null;
    PhaseLR = null;

    // Проверка на допустимость полученных погрешностей
    if(selfTestResult->MaxDiff_ACH                     <= MAX_FFT_DIFF &&
            selfTestResult->MaxDiff_ALG_to_EXP_to_ALG       <= MAX_FFT_DIFF &&
            selfTestResult->MaxDiff_FORWARD_BACKWARD        <= MAX_FFT_DIFF &&
            selfTestResult->MaxDiff_FORWARD_BACKWARD_AntiTW <= MAX_FFT_DIFF &&
            selfTestResult->MaxDiff_PhaseLR                 <= MAX_FFT_DIFF)
    {
        selfTestResult->AllOK = 1;
    }
    else
    {
        selfTestResult->AllOK = 0;
    }

    // Если активирован режим сброса дампа...
    /*if(IsDumpMode)
    {
        // Вердикт по точности выполнения самодиагностики
        DebugHelper.WriteInt(DumpName, "AllOK.int32", selfTestResult->AllOK);

        // Максимальная невязка по расчету заданной АЧХ
        DebugHelper.WriteDouble(DumpName, "MaxDiff_ACH.double", selfTestResult->MaxDiff_ACH);

        // Max. невязка ALG . EXP и обратно
        DebugHelper.WriteDouble(DumpName, "MaxDiff_ALG_to_EXP_to_ALG.double", selfTestResult->MaxDiff_ALG_to_EXP_to_ALG);

        // Max. невязка FORVARD + BACKWARD
        DebugHelper.WriteDouble(DumpName, "MaxDiff_FORWARD_BACKWARD.double", selfTestResult->MaxDiff_FORWARD_BACKWARD);

        //...то же + восст. после TW
        DebugHelper.WriteDouble(DumpName, "MaxDiff_FORWARD_BACKWARD_AntiTW.double", selfTestResult->MaxDiff_FORWARD_BACKWARD_AntiTW);

        // Макс. невязка по расчету разности хода фаз
        DebugHelper.WriteDouble(DumpName, "MaxDiff_PhaseLR.double", selfTestResult->MaxDiff_PhaseLR);

        // Время работы CFFT_Process()
        DebugHelper.WriteDouble(DumpName, "CFFT_Process_time.double", selfTestResult->CFFT_Process_time);

        // Время работы CFFT_Explore()
        DebugHelper.WriteDouble(DumpName, "CFFT_Explore_time.double", selfTestResult->CFFT_Explore_time); 
    }*/

    // Возвращаем результаты тестирования
    return selfTestResult;
}
//-----------------------------------------------------------------------------------------------------------------------------
