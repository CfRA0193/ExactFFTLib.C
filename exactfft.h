#define _USE_MATH_DEFINES
#include <math.h>

#define TRUE 1
#define FALSE 0

typedef char bool;

const double M_PI            = 3.14159265358979323846; //..264338328
const double M_2PI           = 2 * M_PI;
const double FLOAT_MIN       = 3.4E-38;    // Допустимый минимум для операций float
const double MAX_FFT_DIFF    = 1E-7;       // Максимальная погрешность FFT
const double MIN_FRAME_WIDTH = 8;          // Наименьший "рабочий" размер окна FFT
const double MAX_KAISER_BETA = 28;         // Max. Beta (Kaiser window), SL: ~240 dB
const double MAX_PATH        = 256;        // Максимальная длина пути        

const bool DIRECT                 = true;  // Обозначение прямого прохода FFT
const bool REVERSE                = false; // Обозначение обратного прохода FFT
const bool USING_NORM             = true;  // Используется нормализация
const bool NOT_USING_NORM         = false; // Не используется нормализация
const bool USING_TAPER_WINDOW     = true;  // Используется взвешивающее окно
const bool NOT_USING_TAPER_WINDOW = false; // Не используется взвешивающее окно
const bool USING_POLYPHASE        = true;  // Используется полифазное FFT
const bool NOT_USING_POLYPHASE    = false; // Не используется полифазное FFT

typedef struct CFFT_SelfTestResult
{
    //-------------------------------------------------------------------------
    int AllOK; // Результат тестирования точности
    //-------------------------------------------------------------------------
    double MaxDiff_ACH; // Максимальная невязка по расчету заданной АЧХ
    double MaxDiff_ALG_to_EXP_to_ALG; // Max. невязка ALG . EXP и обратно
    double MaxDiff_FORWARD_BACKWARD;  // Max. невязка FORVARD + BACKWARD
    double MaxDiff_FORWARD_BACKWARD_AntiTW; //...то же + восст. после TW
    double MaxDiff_PhaseLR;   // Макс. невязка по расчету разности хода фаз
    double CFFT_Process_time; // Время работы CFFT_Process()
    double CFFT_Explore_time; // Время работы CFFT_Explore()
    //-------------------------------------------------------------------------
} CFFT_SeltTestResult_t;

// Типы взвешивающих окон FFT
// Характеристики окон: PS - "Peak Sidelobe" (наивысший боковой лепесток, дБ)        
enum CosTW
{
    NONE,
    RECTANGULAR_13dbPS,
    HANN_31dbPS,
    HAMMING_43dbPS,
    MAX_ROLLOFF_3_TERM_46dbPS,
    BLACKMAN_58dbPS,
    COMPROMISE_3_TERM_64dbPS,
    EXACT_BLACKMAN_68dbPS,
    MIN_SIDELOBE_3_TERM_71dbPS,
    MAX_ROLLOFF_4_TERM_60dbPS,
    COMPROMISE1_4_TERM_82dbPS,
    COMPROMISE2_4_TERM_93dbPS,
    BLACKMAN_HARRIS_92dbPS,
    NUTTALL_93dbPS,
    BLACKMAN_NUTTALL_98dbPS,
    ROSENFIELD
};
//-----------------------------------------------------------------------------------------------------------------------------
typedef struct CFFT_Object
{
    int N;           // Количество точек FFT
    int NN;          // Кол-во чисел (re + im) FFT
    int NPoly;       // Пересчитанное для полифазного FFT количество точек
    int NNPoly;      // Кол-во чисел (re + im) полифазного FFT
    double Beta;     // Формирующий коэффициент "Beta" окна Кайзера
    CosTW CosTW;     // Тип косинусного взвешивающего окна (если нет - используется окно Кайзера)
    int PolyDiv;     // Делитель "полифазности" FFT ("0" - обычное FFT)
    int WindowStep;  // Шаг окна FFT
    bool IsComplex;  // Используется комплексный вход? (бывает COMPLEX и L+R)
    //-------------------------------------------------------------------------
    int* FFT_P;     // Вектор изменения порядка следования данных перед FFT     // SIZE !!!
    int* FFT_PP;    // Вектор изменения порядка... (для полифазного FFT)        //
    double* FFT_TW; // Взвешивающее окно                                        //
} CFFT_Object_t;
//-----------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Логарифм по произвольному основанию
/// </summary>
/// <param name="arg"> Аргумент логарифма. </param>
/// <param name="logBase"> Основание логарифма. </param>
//-----------------------------------------------------------------------------------------------------------------------------
double LogX(double arg, double logBase);
/// <summary>
/// Приведение значения к ближайшей снизу степени двойки
/// </summary>
/// <param name="arg"> Входной аргумент. </param>
int ToLowerPowerOf2(int arg);
//-----------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Переход от линейной шкалы к dB
/// </summary>
/// <param name="arg"> Входной аргумент (линейная шкала). </param>
/// <param name="zero_db_level"> "Нулевой" уровень dB. </param>
double To_dB(double arg, double zero_db_level);
//-----------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Переход от dB к линейной шкале
/// </summary>
/// <param name="arg"> Входной аргумент (логарифмическая шкала). </param>
/// <param name="zero_db_level"> "Нулевой" уровень dB. </param>
double From_dB(double arg, double zero_db_level);
//-----------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Получение частоты заданного узла FFT
/// </summary>
/// <param name="FFT_Node"> Номер гармоники. </param>
/// <param name="sampFreq"> Частота семплирования. </param>        
/// <param name="N"> Размер окна FFT. </param>
/// <param name="isComplex"> Комплексный режим? </param>
double FreqNode(double FFT_Node, double sampFreq, int N, bool isComplex);
//-----------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Получение узла FFT по заданной частоте
/// </summary>
/// <param name="freqNode"> Заданная частота. </param>
/// <param name="sampFreq"> Частота семплирования. </param>        
/// <param name="N"> Размер окна FFT. </param>
/// <param name="isComplex"> Комплексный режим? </param>
double FFT_Node(double freqNode, double sampFreq, int N, bool isComplex);
//-----------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Нормирование фазы
/// </summary>
/// <param name="phase"> Значение фазы для нормирования. </param>
double PhaseNorm(double phase);
//-----------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Безопасное в смысле значений аргументов вычисление арктангенса
/// </summary>
/// <param name="im"> Мнимая часть комплексного числа. </param>
/// <param name="re"> Действительная часть комплексного числа. </param>
double Safe_atan2(double im, double re);
/// <summary>
/// Заполнение вектора изменения порядка следования данных перед FFT
/// </summary>
/// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
//-----------------------------------------------------------------------------------------------------------------------------
void fill_FFT_P(CFFT_Object* fftObj);
/// <summary>
/// Заполнение вектора изменения порядка следования данных перед FFT
/// (для полифазного FFT)
/// </summary>
/// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
void fill_FFT_PP(CFFT_Object* fftObj);
//-----------------------------------------------------------------------------------------------------------------------------
// Заполнение вектора косинусного взвешивающего окна
void fill_FFT_TW_Cosine(CFFT_Object_t* pfftObj);
//-----------------------------------------------------------------------------------------------------------------------------
// Модифицир. функция Бесселя нулевого порядка первого рода
double BesselI0(double arg);
//-----------------------------------------------------------------------------------------------------------------------------
// Метод заполнения взвешивающего окна (окно Кайзера)
void fill_FFT_TW_Kaiser(CFFT_Object_t* pfftObj);
//-----------------------------------------------------------------------------------------------------------------------------
// "Контролер" объекта FFT
bool CFFT_Inspector(CFFT_Object_t* pfftObj);
//-----------------------------------------------------------------------------------------------------------------------------
void CFFT_Destructor(CFFT_Object* pfftObj);
//-----------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Фабрика объектов FFT
/// </summary>
/// <param name="frameWidth"> Размер кадра. </param>
/// <param name="cosTW"> Тип косинусного взвешивающего окна (если нет - используется окно Кайзера). </param>
/// <param name="beta"> Формирующий коэффициент окна Кайзера. </param>
/// <param name="polyDiv2"> Уровень полифазности как степень двойки. </param>
/// <param name="windowStep"> Шаг окна FFT. </param>
/// <param name="isComplex"> Используется комплексный вход? (бывает COMPLEX и L+R). </param>
CFFT_Object_t* CFFT_Init(int frameWidth, CosTW cosTW, double beta, int polyDiv2,
        int windowStep, bool isComplex);
//-----------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Фабрика объектов FFT
/// </summary>
/// <param name="frameWidth"> Размер кадра. </param>
/// <param name="cosTW"> Тип косинусного взвешивающего окна. </param>
/// <param name="polyDiv2"> Уровень полифазности как степень двойки. </param>
/// <param name="windowStep"> Шаг окна FFT. </param>
/// <param name="isComplex"> Используется комплексный вход? (бывает COMPLEX и L+R). </param>
CFFT_Object* CFFT_Constructor_Cosine(int frameWidth, CosTW cosTW, int polyDiv2, int windowStep, bool isComplex);
//-----------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Фабрика объектов FFT
/// </summary>
/// <param name="frameWidth"> Размер кадра. </param>
/// <param name="beta"> Формирующий коэффициент окна Кайзера. </param>
/// <param name="polyDiv2"> Уровень полифазности как степень двойки. </param>
/// <param name="windowStep"> Шаг окна FFT. </param>
/// <param name="isComplex"> Используется комплексный вход? (бывает COMPLEX и L+R). </param>
CFFT_Object* CFFT_Constructor_Kaiser(int frameWidth, double beta, int polyDiv2, int windowStep, bool isComplex);
//-----------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Основной метод комплексного FFT
/// </summary>
/// <param name="FFT_S"> Вектор входных данных
/// ("левый" и "правый" каналы - чет./нечет.). </param>
/// <param name="FFT_S_Offset"> Смещение данных для анализа во
/// входном векторе FFT_S. </param>
/// <param name="FFT_T"> Выходной вектор коэффициентов. </param>
/// <param name="useTaperWindow"> Использовать взвешивающее окно? </param>
/// <param name="recoverAfterTaperWindow"> Аннигилировать действие
/// взвешивающего окна на обратном проходе? </param>
/// <param name="useNorm"> Использовать нормализацию 1/N? </param>
/// <param name="direction"> Направление преобразования (true - прямое).
/// </param>
/// <param name="usePolyphase"> Использовать полифазное FFT? </param>
/// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
void CFFT_Process(double* FFT_S, int FFT_S_Offset, double* FFT_T, bool useTaperWindow, bool recoverAfterTaperWindow,
        bool useNorm, bool direction, bool usePolyphase, CFFT_Object* pfftObj);
//-----------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Исследование "левого" и "правого" каналов: ("левый" -
/// действительная часть исходных данных, "правый" - мнимая часть)
/// в режиме (L+R)
/// </summary>
/// <param name="FFT_T"> Выходной вектор коэффициентов. </param>
/// <param name="MagL"> Магнитуды "левого" канала. </param>
/// <param name="MagR"> Магнитуды "правого" канала. </param>
/// <param name="ACH"> АЧХ (отношение магнитуды "правого" канала к магнитуде
/// "левого" - как "выход" / "вход"). </param>
/// <param name="ArgL"> Аргумент "левого" канала. </param>
/// <param name="ArgR"> Аргумент "правого" канала. </param>
/// <param name="PhaseLR"> Разность хода фаз каналов ("правый" минус
/// "левый"). </param>
/// <param name="usePolyphase"> Использовать полифазное FFT? </param>
/// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
void CFFT_Explore(double* FFT_T, double* MagL, double* MagR, double* ACH, double* ArgL, double* ArgR,
        double* PhaseLR, bool usePolyphase, CFFT_Object* pfftObj);
//-----------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Исследование результатов комплексного FFT (идентично CFFT из MathCAD)
/// </summary>
/// <param name="FFT_T"> Выходной вектор коэффициентов. </param>
/// <param name="Mag"> Магнитуды. </param>
/// <param name="Arg"> Аргументы. </param>
/// <param name="usePolyphase"> Использовать полифазное FFT? </param>
/// <param name="isMirror"> Зеркальное отображение спектра? </param>
/// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
void CFFT_ComplexExplore(double* FFT_T, double* Mag, double* Arg,
        bool usePolyphase, bool isMirror, CFFT_Object* pfftObj);
//-----------------------------------------------------------------------------------------------------------------------------
/// <summary>
/// Перевод значений массива double в форму dB
/// </summary>
/// <param name="data"> Данные для обработки. </param>
/// <param name="zero_db_level"> Уровень "нулевого" уровня. </param>
/// <param name="fftObj"> Объект FFT, для которого вызывается функция. </param>
void dB_Scale(double* Mag, double zero_db_level, CFFT_Object* pfftObj);
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
CFFT_SelfTestResult SelfTest_S(double* FFT_S, double ACH_Difference, CFFT_Object* pfftObj);
//-----------------------------------------------------------------------------------------------------------------------------
