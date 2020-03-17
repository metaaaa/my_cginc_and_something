Shader "Unlit/raymarching_test"
{
Properties
{
_MainTex ("Texture", 2D) = "white" {}
_Value1 ("Value1", Float) = 1.0
_Value2 ("Value2", Float) = 1.0
_Value3 ("Value3", Float) = 1.0
_Value4 ("Value4", Float) = 1.0
_Range1 ("Range1", Range (-50, 50)) = 0
_Color1 ("Color1", Color) = (1.0,1.0,1.0,1.0)
_Vector1 ("Vector1", Vector) = (0.0,0.0,0.0,0.0)
}
SubShader
{
Tags { "RenderType"="Opaque" }
LOD 100
Cull Front

Pass
{
CGPROGRAM
#pragma vertex vert
#pragma fragment frag

#include "UnityCG.cginc"

#include "meta_function.cginc"
#include "meta_DistanceFunctions.cginc"


struct appdata
{
    float4 vertex : POSITION;
};

struct v2f
{
    float4 vertex : SV_POSITION;
    float3 worldPos : TEXCOORD1;
};

struct f_out
{
    fixed4 color : SV_Target;
    float depth : SV_Depth;
};

sampler2D _MainTex;
float4 _MainTex_ST;
float _Value1;
float _Value2;
float _Value3;
float _Value4;
float _Range1;
fixed4 _Color1;
float4 _Vector1;

float distanceFunc(float3 p)
{
    //p.xy=rot(p.xy,_Time.y);
    //p.xz=rot(p.xz,_Time.y);
    //p.yz=rot(p.yz,_Time.y);
    p.xy=rot(p.xy,-PI/2.0);
    //p.xy=rot(p.xy,p.z*_Value1);
    
    return deQuaternionMandelbrot(float4(p - float3(0.5, 0.0, 0.0), 0.0));
}

            
// 法線を得る
float3 getNormal(float3 p){
    const float d = 0.0001;
    return normalize(float3(
        distanceFunc(p + float3(  d, 0.0, 0.0)) - distanceFunc(p + float3( -d, 0.0, 0.0)),
        distanceFunc(p + float3(0.0,   d, 0.0)) - distanceFunc(p + float3(0.0,  -d, 0.0)),
        distanceFunc(p + float3(0.0, 0.0,   d)) - distanceFunc(p + float3(0.0, 0.0,  -d))
    ));
}

// 深度を得る
float getDepth(float3 pos)
{
    float4 sPos = UnityObjectToClipPos(float4(pos, 1.0));
#if defined(SHADER_TARGET_GLSL)
    return (sPos.z / sPos.w) * 0.5 + 0.5;
#else 
    return sPos.z / sPos.w;
#endif 
}

v2f vert (appdata v)
{
    v.vertex.xyz*=4;
    v2f o;
    o.vertex = UnityObjectToClipPos(v.vertex);
    o.worldPos = v.vertex.xyz;
    return o;
}

f_out frag (v2f i)
{
    f_out o;

    float3 rayOrg = mul(unity_WorldToObject, float4(_WorldSpaceCameraPos, 1)).xyz;

    float3 rayDir = normalize(i.worldPos - rayOrg);

    float3 rayPos = rayOrg;
    
    float distance = 0.0;
    for(int i = 0; i < 32; i++)
    {
        distance = distanceFunc(rayPos);
        rayPos += distance * rayDir;
    }
    
    if(abs(distance) < 0.01)
    {
        float3 normal = getNormal(rayPos);
        normal = normalize(UnityObjectToWorldNormal(normal));
        o.color.xyz = 0.5 + 0.5 * normal;
        o.color.w = 1.0;
    }
    else
    {
        discard;
    }

    o.depth = getDepth(rayPos);

    return o;
}
ENDCG
}
}
}
