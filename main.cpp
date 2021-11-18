#include <iostream>
#include <gmp.h>
#include <cstring>
#include <sstream>

namespace ElGamal {
#define MSGLEN 100
    unsigned char s[100000];
    mpz_t q, p, alpha, d, beta, k;
    mpz_t p_1,t1,t2;
    mpz_t c1,c2[MSGLEN];
    gmp_randstate_t state;
    size_t key_size, nc;
    char f[16]=  {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};
    char f1[256];

    //下面四个函数调试用
    template<typename T> void hexDump(const std::string& info, T x){
        int sx = (sizeof x)*8-4;
        std::cout<<info<<' ';
        for (int i=sx; i>=0; i-=4){
            putchar(f[(x>>i)&0xf]);
        }
        putchar('\n');
    }
    template<typename T> void hexDump(const std::string& info, T *x, size_t length){
        std::cout<<info<<' ';
        if (length < 1) return;
        int sx = (sizeof x[0])*8-4;
        for  (int i=0; i<length; i++){
            for (int j=sx; j>=0; j-=4){
                putchar(f[(x[i]>>j)&0xf]);
            }
            putchar(' ');
        }
        putchar('\n');
    }
    template<typename T> void binaryDump(const std::string& info,T x){
        int bits = (sizeof x) * 8;
        std::cout<<info<<' ';
        for (int i= bits - 1; i >= 0; i--){
            putchar(((x>>i)&1)?'1':'0');
            if ((i & 7) == 0) putchar(' ');
        }
        putchar('\n');
    }
    template<typename T> void binaryDump(const std::string& info, T *x, size_t length){
        std::cout<<info<<' ';
        if (length < 1) return;
        int bits = (sizeof x[0]) * 8;
        for (size_t i = 0; i<length; i++){
            for (int j=bits-1; j>=0; j--){
                putchar(((x[i]>>j)&1)?'1':'0');
                if ((j & 7) == 0) putchar(' ');
            }
            putchar(' ');
        }
        putchar('\n');
    }

    /**
     * 打印菜单，用于用户选择密钥长度
     */
    void menu(){
        printf(""
				" _____   _    ____                               _ \n"
        		"| ____| | |  / ___|   __ _   _ __ ___     __ _  | |\n"
        		"|  _|   | | | |  _   / _` | | '_ ` _ \\   / _` | | |\n"
        		"| |___  | | | |_| | | (_| | | | | | | | | (_| | | |\n"
        		"|_____| |_|  \\____|  \\__,_| |_| |_| |_|  \\__,_| |_|"
        );
        printf(""
               "\n============================================="
               "\n|  Please choose size of encryption size    |"
               "\n|         1.  64 bits                       |"
               "\n|         2.  128 bits                      |"
               "\n|         3.  256 bits                      |"
               "\n|         4.  512 bits                      |"
               "\n|         5.  1024 bits                     |"
               "\n|         6.  2048 bits                     |"
               "\n============================================="
               "");
        int n = -1;
        int choice;
        while(n == -1){
            printf("\n Please input your choice:");
            scanf("%d",&choice);
            getchar();
            switch(choice) {
                case 1:n=64;  break;
                case 2:n=128; break;
                case 3:n=256; break;
                case 4:n=512; break;
                case 5:n=1014;break;
                case 6:n=2048;break;
                default :
                    n  = -1;
                    printf("Sorry this is not an available option\n\n");
            }
        }
        key_size = n;
    }

    /**
     * 初始工作
     * 初始化大整数对象和符号逆变换
     */
    void init(){
        mpz_init(q);mpz_init(p);mpz_init(alpha);
        mpz_init(d);mpz_init(beta);mpz_init(k);
        mpz_init(p);mpz_init(p_1);mpz_init(t1);
        mpz_init(t2);mpz_init(c1);
        //gmp_randinit_default(state);
        //用于随机数生成
        gmp_randinit_mt(state);
        //根据系统时间设置random的seed
        gmp_randseed_ui(state,(unsigned int)(time(NULL)));
        for (int i=0; i<MSGLEN; i++){
            mpz_init2(c2[i], key_size);
        }
        // 字符到数值逆转换
        f1['1'] = 1;  f1['2'] = 2;
        f1['3'] = 3;  f1['4'] = 4;
        f1['5'] = 5;  f1['6'] = 6;
        f1['7'] = 7;  f1['8'] = 8;
        f1['9'] = 9;
        f1['A'] = 10; f1['a'] = 10;
        f1['B'] = 11; f1['b'] = 11;
        f1['C'] = 12; f1['c'] = 12;
        f1['D'] = 13; f1['d'] = 13;
        f1['E'] = 14; f1['e'] = 14;
        f1['F'] = 15; f1['f'] = 15;
    }

    /**
     * 生成密匙
     * 主要是p、alpha、d、beta
     */
    void keyGen(){
        while(true) {
            //生成key_size位的随机数
            mpz_urandomb(q, state, key_size);
            mpz_nextprime(q, q);//生成素数
            // p = 2*q + 1，构造一个安全素数p
            // p-1因子只有1，2，q，p-1四个
            mpz_mul_ui(p, q, 2);
            mpz_add_ui(p, p, 1);
            // 使用GMP内置的Miller-Rabin测试p是否素数
            // 当key_size过大时，程序将会花费大量时间判断
            // 生成的p是否素数
            // 简易使用 key_size < 1024的情况
            if (mpz_probab_prime_p(p,15)>0) break;
        }
        gmp_printf("p :     %Zd\n", p);
        // p_1 = p-1,用于alpha、d生成
        mpz_sub_ui(p_1, p, 1);
        while(true){//寻找p的本原元 alpha
            mpz_urandomm(alpha, state, p_1);//生[0,p-2]的随机数
            if (mpz_cmp_ui(alpha, 1) == 0){// alpha = 1，重新选择
                continue;
            }
            mpz_powm_ui(t1, alpha, 2, p);//t1 = alpha^2 % p
            mpz_powm(t2, alpha, q, p);   //t2 = alpha^q % p
            if (mpz_cmp_ui(t1, 1) != 0 && mpz_cmp_ui(t2, 1) != 0){
                // g^2 % p != 1 and g^q % p != 1
                break;
            }
        }
        gmp_printf("alpha : %Zd\n",alpha);
        while(true){//生成随机整数d
            mpz_urandomm(d, state, p_1);
            if (mpz_cmp_ui(d, 0)!=0){//确保d!=0
                break;
            }
        }
        gmp_printf("d :     %Zd\n",d);
        // beta = alpha^d % p
        mpz_powm(beta, alpha, d, p);
        gmp_printf("beta :  %Zd\n",beta);
    }

    /**
     * 将src按位拓展，只能分成2、4、8份
     * @param src 源字符串
     * @param len_src 源字符串长度
     * @param mul 扩长倍数
     * @return 返回扩长后的字符指针（函数内部自行malloc）
     */
    unsigned char * strExpand(unsigned char* src, size_t len_src, int mul=2){
        static int b[]=  {0,0,4,  0,6,0,0,0,7};
        static int det[]={0,0,4,  0,2,0,0,0,1};
        static int msk[]={0,0,0xF,0,3,0,0,0,1};
        if (mul != 2 && mul != 4 && mul != 8) return nullptr;
        auto *des = (unsigned char*)malloc(len_src*mul*sizeof(char));
        for (int i=0; i<len_src; i++){
            for (int j=b[mul],kk=0; j>=0; j-=det[mul],kk++){
                des[i*mul+kk] = f[((src[i])>>j)&msk[mul]];
            }
        }
        des[len_src*mul] = 0;
        return des;
    }

    /**
     * strExpand的逆变换
     * @param src 源字符串
     * @param len_src 源字符串长度
     * @param mul 扩长倍数
     * @return 返回扩长前的字符指针（函数内部自行malloc）
     */
    unsigned char * strNarrow(unsigned char *src, size_t len_src, int mul=2){
        size_t malloc_size = (((len_src+mul-1)/mul)+1)*sizeof(char);
        auto des = (unsigned char *)malloc(malloc_size);
        int len = 0;
        des[0] = 0;
        for (int i=0; i<len_src; i++){
            if (i!=0 && i%2 == 0){
                len++;
                des[len]=0;
            }
            des[len] = ((des[len]<<4)|f1[src[i]]);
        }
        des[++len] = 0;
        return des;
    }

    /**
     * 加密入口
     * @param str 要加密的字符串
     * @param len 字符串长度
     */
    void encryption(char *str, size_t len){
        while(true){//生成 [1,p-1] 的随机数 k
            mpz_urandomm(k, state, p_1);
            if (mpz_cmp_ui(k, 0)!=0){//确保k大于0
                break;
            }
        }
        // c1 = alpha^k % p
        mpz_powm(c1, alpha, k, p);
        gmp_printf("c1 :    %Zd\n",c1);
        // t1 = beta^k % p
        mpz_powm(t1, beta, k, p);

        /* 出于性能和安全性考虑，需将特定长度字符串压缩为一个大整数
         * 1）在ElGamal加密中，明文要求小于p
         * 2）通过将字符串4个比特位化为一个16进制数，才能方便使用
         *  函数 mpz_set_str 将数存入 gmp 大整数对象。这里近似
         *  于将明文“分块”，每一块的大小是((key_size/8)-1)*2
         *  （乘以2是因为明文通过strExpand扩展为两倍大小）
         */
        auto tstr = strExpand((unsigned char*)str,len);
        int blocks = ((key_size/8)-1)*2;
        nc = 0;
        unsigned char ch;
        printf("====================================\n"
               "|convert char array to big integers|\n"
               "====================================\n");
        for (int i=0; i<len*2; i+=blocks){
            if (i+blocks < len*2){
                ch = tstr[i+blocks];
                tstr[i+blocks] = 0;
            }
            //[tstr+i,tsrt+blocks)区间构成一个大整数
            mpz_set_str(c2[nc], reinterpret_cast<const char *>(tstr + i), 16);
            mpz_get_str((char *)s, 16, c2[nc]);
            gmp_printf("m[%d] :  %Zd\n",nc, c2[nc]);

            //c2 = c2*t1 % p = m * beta^k % p
            mpz_mul(c2[nc], c2[nc], t1);
            mpz_mod(c2[nc], c2[nc], p);
            nc++;
            if (i+blocks < len*2){
                tstr[i+blocks] = ch;
            }
        }
        free(tstr);
        printf("====================================\n"
               "|do encryption                     |\n"
               "====================================\n");
        for (int i=0; i<nc; i++){
            gmp_printf("c2[%d] : %Zd\n",i, c2[i]);
        }
    }

    /**
     * 解密函数
     */
    void decryption(){
        //存储最终的流结果
        std::stringstream sstream;
        // c1 = (c1^d)^-1 (mod p)
        mpz_powm(c1,c1,d,p);
        mpz_invert(c1,c1,p);//求逆元
        printf("====================================\n"
               "|do decryption                     |\n"
               "====================================\n");
        for (int i=0; i<nc; i++){
            //m = c2 * c1 % p = c2 * (c1^d)^-1 % p
            mpz_mul(c2[i],c2[i],c1);
            mpz_mod(c2[i],c2[i],p);
            gmp_printf("m[%d] :  %Zd\n",i,c2[i]);
            mpz_get_str((char*)s,16,c2[i]);
            //字符串压缩：两位变一位（16进制）得到明文消息
            auto ss = strNarrow(s,strlen((char*)s));
            sstream<<ss;
        }
        std::cout<<sstream.str()<<std::endl;
    }

    /**
     * ElGamal入口
     * @param src 要加密的字符串
     * @param len 要加密的长度
     * 程序自动生成随机数、密钥等
     */
    void main(char *src, int len){
        menu();
        init();
        keyGen();
        encryption(src, len);
        decryption();
    }
}

int main() {
    char ss[] = "zhouzhiiiiiiiidiiiiiiziiiiioiiiiiiiipiiiiiiiiiiiiiiciiiik";
    ElGamal::main(ss, strlen(ss));
    return 0;
}
