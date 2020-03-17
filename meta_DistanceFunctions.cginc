float sdSphere( float3 p, float s )
{
  return length(p)-s;
}

float sdBox( float3 p, float3 b )
{
  float3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float sdRoundBox( float3 p, float3 b, float r )
{
  float3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0) - r;
}

float sdPlane( float3 p, float4 n )
{
  // n must be normalized
  return dot(p,n.xyz) + n.w;
}

float sdHexPrism( float3 p, float2 h )
{
  const float3 k = float3(-0.8660254, 0.5, 0.57735);
  p = abs(p);
  p.xy -= 2.0*min(dot(k.xy, p.xy), 0.0)*k.xy;
  float2 d = float2(
       length(p.xy-float2(clamp(p.x,-k.z*h.x,k.z*h.x), h.x))*sign(p.y-h.x),
       p.z-h.y );
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdTriPrism( float3 p, float2 h )
{
    const float k = sqrt(3.0);
    h.x *= 0.5*k;
    p.xy /= h.x;
    p.x = abs(p.x) - 1.0;
    p.y = p.y + 1.0/k;
    if( p.x+k*p.y>0.0 ) p.xy=float2(p.x-k*p.y,-k*p.x-p.y)/2.0;
    p.x -= clamp( p.x, -2.0, 0.0 );
    float d1 = length(p.xy)*sign(-p.y)*h.x;
    float d2 = abs(p.z)-h.y;
    return length(max(float2(d1,d2),0.0)) + min(max(d1,d2), 0.);
}

float sdCapsule( float3 p, float3 a, float3 b, float r )
{
  float3 pa = p - a, ba = b - a;
  float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
  return length( pa - ba*h ) - r;
}

float sdVerticalCapsule( float3 p, float h, float r )
{
  p.y -= clamp( p.y, 0.0, h );
  return length( p ) - r;
}

float sdCylinder( float3 p, float3 c )
{
  return length(p.xz-c.xy)-c.z;
}

float sdCappedCylinder( float3 p, float h, float r )
{
  float2 d = abs(float2(length(p.xz),p.y)) - float2(h,r);
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdCone( float3 p, float2 c )
{
  // c is the sin/cos of the angle
  float q = length(p.xy);
  return dot(c,float2(q,p.z));
}

float sdTorus( float3 p, float2 t )
{
  float2 q = float2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}

float sdLink( float3 p, float le, float r1, float r2 )
{
  float3 q = float3( p.x, max(abs(p.y)-le,0.0), p.z );
  return length(float2(length(q.xy)-r1,q.z)) - r2;
}

float sdOctahedron( float3 p, float s)
{
  p = abs(p);
  return (p.x+p.y+p.z-s)*0.57735027;
}

float sdPyramid( float3 p, float h)
{
  float m2 = h*h + 0.25;
    
  p.xz = abs(p.xz);
  p.xz = (p.z>p.x) ? p.zx : p.xz;
  p.xz -= 0.5;

  float3 q = float3( p.z, h*p.y - 0.5*p.x, h*p.x + 0.5*p.y);
   
  float s = max(-q.x,0.0);
  float t = clamp( (q.y-0.5*p.z)/(m2+0.25), 0.0, 1.0 );
    
  float a = m2*(q.x+s)*(q.x+s) + q.y*q.y;
  float b = m2*(q.x+0.5*t)*(q.x+0.5*t) + (q.y-m2*t)*(q.y-m2*t);
    
  float d2 = min(q.y,-q.x*m2-q.y*0.5) > 0.0 ? 0.0 : min(a,b);
    
  return sqrt( (d2+q.z*q.z)/m2 ) * sign(max(q.z,-p.y));
}

//メンガースポンジ
//return deMengerSponge2(p, float3(_Value,_Value,_Value), _Value2);//0.44,3.1
float deMengerSponge2(float3 p, float3 offset, float scale) 
{
	float4 z = float4(p, 1.0);
	for (int i = 0; i < 5; i++) {
		z = abs(z);
		if (z.x < z.y) z.xy = z.yx;
		if (z.x < z.z) z.xz = z.zx;
		if (z.y < z.z) z.yz = z.zy;
		z *= scale;
		z.xyz -= offset * (scale - 1.0);
		if (z.z < -0.5 * offset.z * (scale - 1.0))
			z.z += offset.z * (scale - 1.0);
		}
	return (length(max(abs(z.xyz) - float3(1.0, 1.0, 1.0), 0.0))) / z.w;
}

//再帰的四面体
//return deRecursiveTetrahedron(p, float3(_Value,_Value,_Value), _Value2);
float deRecursiveTetrahedron(float3 p, float3 offset, float scale) {
    float4 z = float4(p, 1.0);
    for (int i = 0; i < 5; i++) {
        if (z.x + z.y < 0.0) z.xy = -z.yx;
        if (z.x + z.z < 0.0) z.xz = -z.zx;
        if (z.y + z.z < 0.0) z.zy = -z.yz;
        z *= scale;
        z.xyz -= offset * (scale - 1.0);
    }
    return (length(z.xyz) - 1.5) / z.w;
}


float4 qmul(float4 a, float4 b) {
    return float4(
        a.x * b.x - a.y * b.y - a.z * b.z - a.w * b.w,
        a.x * b.y + a.y * b.x - a.z * b.w + a.w * b.z,
        a.x * b.z + a.y * b.w + a.z * b.x - a.w * b.y,
        a.x * b.w - a.y * b.z + a.z * b.y + a.w * b.x
    );
}
//なんかぐるぐるしてるやつ
//return deQuaternionJuliaSet(float4(p, 0.0), float4(-1.0, 0.2, 0.0, 0.0));
float deQuaternionJuliaSet(float4 p, float4 c) {
    float4 z = p;
    float4 dz = float4(1.0, 0.0, 0.0, 0.0);
    float4 pz, pdz;
    float r = 0.0, dr = 1.0;
    for (int i = 0; i < 16; i++) {
        pz = z;
        z = qmul(pz, pz) + c;
        pdz = dz;
        dz = 2.0 * qmul(pz, pdz);
        r = length(z);
        dr = length(dz);
        if (r > 4.0) break;
    }
    return 0.5 * log(r) * r / dr;
}

//噴水みたいなやつ
//return deQuaternionMandelbrot(float4(p - float3(0.5, 0.0, 0.0), 0.0));
float deQuaternionMandelbrot(float4 p) {
    float4 z = float4(0.0,0.,0.,0.);
    float4 dz = float4(0.0,0.,0.,0.);
    float4 pz, pdz;
    float r, dr;
    for (int i = 0; i < 8; i++) {
        pz = z;
        z = qmul(pz, pz) + p;
        pdz = dz;
        dz = 2.0 * qmul(pz, pdz) + 1.0;
        r = length(z);
        dr = length(dz);
        if (r > 8.0)  break;

    }
    return 0.5 * log(r) * r / dr;
}
//六角形タイリング
//return max(deHexTiling(p.zx, 1.0, 0.9), p.y);
float sdHex(float2 p, float h) {
    float3 k = float3(-0.8660254, 0.57735, 0.5);
    p = abs(p);
    p -= 2.0 * min(dot(k.xz, p), 0.0) * k.xz;
    return length(p - float2(clamp(p.x, -k.y * h, k.y * h), h)) * sign(p.y - h);
}

// SQRT3 = sqrt(3.0)
#define SQRT3 1.73205080757
// 0.0 <= scale <= 1.0
float deHexTiling(float2 p, float radius, float scale) {
    float2 rep = float2(2.0 * SQRT3, 2.0) * radius;
    float2 p1 = fmod(p, rep) - rep * 0.5;
    float2 p2 = fmod(p + 0.5 * rep, rep) - rep * 0.5;
    return min(
        sdHex(p1.xy, scale * radius),
        sdHex(p2.xy, scale * radius)
    );
}

//マンデルバルブ
//return deMandelbulb(p, 8.0);
float deMandelbulb(float3 p, float power) {
    float3 z = p;
    float dr = 1.0;
    float r;
    for (int i = 0; i < 8; i++) {
        r = length(z);
        if (r > 10.0) break;
        float theta = acos(z.y / r);
        float phi = atan2(z.z, z.x);
        dr = pow(r, power - 1.0) * power * dr + 1.0;

        float zr = pow(r, power);
        theta = theta * power;
        phi = phi * power;

        z = zr * float3(sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi));
        z += p;
    }
    return 0.5 * log(r) * r / dr;
}

//マンデルボックス
//return deMandelbox(p, 2.0, 0.5, 1.0);
float3 boxFold(float3 z, float dz) {
    return clamp(z, 1.0, 1.0) * 2.0 - z;
}

void sphereFold(inout float3 z, inout float dz, float minRadius, float fixedRadius) {
    float m2 = minRadius * minRadius;
    float f2 = fixedRadius * fixedRadius;
    float r2 = dot(z, z);
    if (r2 < m2) {
        float temp = (f2 / m2);
        z *= temp;
        dz *= temp;
    } else if (r2 < f2) {
        float temp = (f2 / r2);
        z *= temp;
        dz *= temp;
    }
}

// ref: http://blog.hvidtfeldts.net/index.php/2011/11/distance-estimated-3d-fractals-vi-the-mandelbox/
float deMandelbox(float3 p, float scale, float minRadius, float fixedRadius) {
    float3 z = p;
    float dr = 1.0;
    for (int i = 0; i < 12; i++) {
        z = boxFold(z, dr);
        sphereFold(z, dr, minRadius, fixedRadius);
        z = scale * z + p;
        dr = dr * abs(scale) + 1.0;
    }
    float r = length(z);
    return r / abs(dr);
}


//return dePseudoKleinian(p);
float dePseudoKleinian(float3 p) {
    float3 csize = float3(0.90756, 0.92436, 0.90756);
    float size = 1.0;
    float3 c = float3(0.0,0.0,0.0);
    float defactor = 1.0;
    float3 offset = float3(0.0,0.0,0.0);
    float3 ap = p + 1.0;
    for (int i = 0; i < 8; i++) {
        ap = p;
        p = 2.0 * clamp(p, -csize, csize) - p;
        float r2 = dot(p, p);
        float k = max(size / r2, 1.0);
        p *= k;
        defactor *= k;
        p += c;
    }
    float r = abs(0.5 * abs(p.y - offset.y) / defactor);
    return r;
}
