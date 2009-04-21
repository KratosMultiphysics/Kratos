#ifndef COLOR_H
#define COLOR_H

static float color_alpha_default=.2f;

// glMaterial no acepta float

struct color{ // default publico
  union{
    float rgba[4];
    struct{float r,g,b,a;};
  };

  // constructores
  color()
    {r=g=b=0; a=color_alpha_default;}
  color(float ri, float gi, float bi, float ai=color_alpha_default) 
    {r=ri;g=gi;b=bi;a=ai;}
  color(const color& c) 
    {r=c.r;g=c.g;b=c.b;a=c.a;}
  color(const float c[4]) 
    {r=c[0];g=c[1];b=c[2];a=c[3];}
  color(const int c[4]) 
    {r=c[0]/255.0f;g=c[1]/255.0f;b=c[2]/255.0f;a=c[3]/255.0f;}

  // destructor
  ~color() {}

  // cordenadas
  float &operator[](int i)       {return rgba[i];}
  float  operator[](int i) const {return rgba[i];}

  //referencia
  operator float*() {return rgba;}
  operator const float*() const {return rgba;}

  // copia
  color& operator=(const color &c) 
    {r=c.r;g=c.g;b=c.b;a=c.a;return *this;}
  color& operator=(const float c[4]) 
    {r=c[0];g=c[1];b=c[2];a=c[3];return *this;}


  color& aclara(float f)
    {r+=f*(1-r);g+=f*(1-g);b+=f*(1-b);return *this;}
  friend const color aclara(const color &c,float f)
    {return color(c.r+f*(1-c.r),c.g+f*(1-c.g),c.b+f*(1-c.b),c.a);}

  color& oscurece(float f)
    {r*=(1-f);g*=(1-f);b*=(1-f);return *this;}
  friend const color oscurece(const color &c,float f)
    {return color(c.r*(1-f),c.g*(1-f),c.b*(1-f),c.a);}
};

// constantes globales
static const color NEGRO(0.f,0.f,0.f,1.f);
static const color BLANCO(1.f,1.f,1.f,1.f);
static const size_t SZcolor=sizeof(color);

#endif
