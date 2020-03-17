#define PI acos(-1.0)

float2 rot(float2 p, float r)
{
    float c = cos(r);
    float s = sin(r);
    return mul(p, float2x2(c, -s, s, c));
}

float3 rand3D(float3 p)
{
    p = float3( dot(p,float3(127.1, 311.7, 74.7)),
    dot(p,float3(269.5, 183.3,246.1)),
    dot(p,float3(113.5,271.9,124.6)));
    return frac(sin(p)*43758.5453123);
}

float rand(float2 co) 
{
    return frac(sin(dot(co.xy, float2(12.9898, 78.233))) * 43758.5453);
}

float rand(float3 co)
{
    return frac(sin(dot(co.xyz, float3(12.9898, 78.233, 56.787))) * 43758.5453);
}

float noise(float3 pos)
{
    float3 ip = floor(pos);
    float3 fp = smoothstep(0, 1, frac(pos));
    float4 a = float4(
        rand(ip + float3(0, 0, 0)),
        rand(ip + float3(1, 0, 0)),
        rand(ip + float3(0, 1, 0)),
        rand(ip + float3(1, 1, 0)));
    float4 b = float4(
        rand(ip + float3(0, 0, 1)),
        rand(ip + float3(1, 0, 1)),
        rand(ip + float3(0, 1, 1)),
        rand(ip + float3(1, 1, 1)));
    a = lerp(a, b, fp.z);
    a.xy = lerp(a.xy, a.zw, fp.y);
    return lerp(a.x, a.y, fp.x);
}

float perlin(float3 pos) {
    return  (noise(pos) * 32 +
            noise(pos * 2 ) * 16 +
            noise(pos * 4) * 8 +
            noise(pos * 8) * 4 +
            noise(pos * 16) * 2 +
            noise(pos * 32) ) / 63;
}



float3 RGB2HSV(float3 c)
{
        float4 K = float4(0.0, -1.0 / 3.0, 2.0 / 3.0, -1.0);
	float4 p = lerp(float4(c.bg, K.wz), float4(c.gb, K.xy), step(c.b, c.g));
	float4 q = lerp(float4(p.xyw, c.r), float4(c.r, p.yzx), step(p.x, c.r));

	float d = q.x - min(q.w, q.y);
	float e = 1.0e-10;
	return float3(abs(q.z + (q.w - q.y) / (6.0 * d + e)), d / (q.x + e), q.x);
}

float3 HSV2RGB(float3 c)
{
      float4 K = float4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
      float3 p = abs(frac(c.xxx + K.xyz) * 6.0 - K.www);
      return c.z * lerp(K.xxx, saturate(p - K.xxx), c.y);
}

//更新頻度を下げる奴、FPSを下げる？
float posterize(float f,float c)
{
		float pstr=floor(f*c)/(c-1);
		return pstr;
}
