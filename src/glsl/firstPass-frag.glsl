/* 
	Code written by Joseph Klinger and Tabatha Hickman
	University of Pennsylvania
	CIS 700 Instructor: Rachel Hwang
	May 2017
*/

#define SPHERE_TRACING true
#define T_MAX 10.0

varying vec2 f_uv;

uniform float u_time;
uniform vec2 u_resolution;
uniform float u_fovy;
uniform float u_aspect;

uniform mat4 u_cwMat;
uniform mat4 u_ccwMat;
uniform mat4 u_northMat;
uniform mat4 u_southMat;
uniform mat4 u_westMat;
uniform mat4 u_eastMat;
uniform mat4 u_rotateX1;
uniform mat4 u_rotateY1;
uniform mat4 u_rotateZ1;
uniform mat4 u_rotateX2;
uniform mat4 u_rotateY2;
uniform mat4 u_rotateZ2;


vec4 resColor;

/***** Geometry SDF Functions
http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
							  											*****/

float SDF_Sphere( vec3 pos, float radius ) {
	return length(pos) - radius;
}

// Taken from IQ's realtime ShaderToy implementation here: https://www.shadertoy.com/view/ltfSWn
float SDF_Mandelbulb( vec3 p , float manPower)
{
	vec3 w = p;
    float m = dot(w,w);

    vec4 trap = vec4(abs(w),m);
    float dz = 1.0;
    
    
    for( int i = 0; i < 4; i++ )
    {
#if 1
        float m2 = m*m;
        float m4 = m2*m2;
        dz = manPower*sqrt(m4*m2*m)*dz + 1.0;

        float x = w.x; float x2 = x*x; float x4 = x2*x2;
        float y = w.y; float y2 = y*y; float y4 = y2*y2;
        float z = w.z; float z2 = z*z; float z4 = z2*z2;

        float k3 = x2 + z2;
        float k2 = inversesqrt( k3*k3*k3*k3*k3*k3*k3 );
        float k1 = x4 + y4 + z4 - 6.0*y2*z2 - 6.0*x2*y2 + 2.0*z2*x2;
        float k4 = x2 - y2 + z2;

        w.x = p.x +  64.0*x*y*z*(x2-z2)*k4*(x4-6.0*x2*z2+z4)*k1*k2;
        w.y = p.y + -16.0*y2*k3*k4*k4 + k1*k1;
        w.z = p.z +  -8.0*y*k4*(x4*x4 - 28.0*x4*x2*z2 + 70.0*x4*z4 - 28.0*x2*z2*z4 + z4*z4)*k1*k2;
#else
        dz = 8.0*pow(m,3.5)*dz + 1.0;
        
        float r = length(w);
        float b = 8.0*acos( clamp(w.y/r, -1.0, 1.0));
        float a = 8.0*atan( w.x, w.z );
        w = p + pow(r,8.0) * vec3( sin(b)*sin(a), cos(b), sin(b)*cos(a) );
#endif        
        
        trap = min( trap, vec4(abs(w),m) );

        m = dot(w,w);
        if( m > 4.0 )
            break;
    }
    trap.x = m;
    resColor = trap;

    return 0.25 * log(m) * sqrt(m) / dz;
}

// SDF Operators:
float intersection(float d1, float d2) {
    return max(d1,d2);
}

float subtraction( float d1, float d2 ) {
    return max(-d1,d2);
}

float un(float d1, float d2) {
    return min(d1,d2);
}

// Returns transformed point based on rotation and translation matrix of shape
vec3 transform( vec3 point, mat4 trans ) {
	//columns of the rotation matrix transpose
	vec3 col1 = vec3(trans[0][0], trans[1][0], trans[2][0]);
	vec3 col2 = vec3(trans[0][1], trans[1][1], trans[2][1]);
	vec3 col3 = vec3(trans[0][2], trans[1][2], trans[2][2]);

	mat3 rotTranspose = mat3(col1, col2, col3);

	vec3 col4 = -1.0*rotTranspose*vec3(trans[3]);

	mat4 newTrans = mat4(vec4(col1, 0.0), vec4(col2, 0.0), vec4(col3, 0.0), vec4(col4, 1.0));

	return vec3(newTrans * vec4(point, 1.0));
}

float mod(int num1, int num2)
{
	int div = num1/num2;
	return float(num1 - div*num2);
}

// Decide which scene to display (time - dependent)
int sceneNum()
{
	float x = u_time;
	float cycle = 124.0;
	float fps = 6.0;
	float t = mod(x, cycle);
	if(t <= 42.0) { return 2; }
	else if(t <= 80.0) { return 1; }
	else { return 3; }
}

// Estimate the distance to the objects in the scene depending on the sceneNum() (see above)
float sceneMap( vec3 pos ) {
	float t = u_time/4.0;
	int sceneNumber = sceneNum();

	if(sceneNumber == 1)
	{
		//SCENE 01------------------------------------------------------------
		float dist1;
		vec3 newPos1 = transform(transform(pos + vec3(sin(t)*3.25, sin(t)*2.0, cos(t)*3.25), u_cwMat), u_northMat);	
		float bb1 = SDF_Sphere(newPos1, 1.1);
		if(bb1 < .015)
		{
			float power = 10.0;
			dist1 = SDF_Mandelbulb(newPos1, power);
		}	
		else
		{
			dist1 = bb1;
		}

		float dist2;
		vec3 newPos2 = transform(transform(pos + vec3(sin(t + 30.0)*3.25, cos(t + 8.0)*2.0, sin(t)*-1.5), u_ccwMat), u_eastMat);	
		float bb2 = SDF_Sphere(newPos2, 1.1);
		if(bb2 < .015)
		{		
			float power = 10.0;
			dist2 = SDF_Mandelbulb(newPos2, power);
		}
		else
		{
			dist2 = bb2;
		}

		float dist3;
		vec3 newPos3 = transform(transform(pos + vec3(cos(t+6.0)*3.25, -1.0*sin(t) + -0.5*cos(1.0), sin(t+12.0)*3.25), u_cwMat), u_westMat);	
		float bb3 = SDF_Sphere(newPos3, 1.1);
		if(bb3 < .015)
		{
		
			float power = 10.0;
			dist3 = SDF_Mandelbulb(newPos3, power);
		}
		else
		{
			dist3 = bb3;
		}

		return un(dist1, un(dist2, dist3));
	}
	else if(sceneNumber == 2)
	{
		//SCENE 02------------------------------------------------------------
		float dist4;
		float displace = pow(mod(u_time, 124.0)/(42.0), log(0.2) / log(0.5)) * 3.0; 
		vec3 newPos4 = transform(transform(transform(pos + vec3(displace, 0, displace), u_rotateY1), u_rotateZ1), u_rotateX1);	
		float bb4 = SDF_Sphere(newPos4, 1.1);
		if(bb4 < .015)
		{	
			float power = 10.0;
			dist4 = SDF_Mandelbulb(newPos4, power);
		}
		else
		{
			dist4 = bb4;
		}

		return dist4;
	}
	else 
	{
		//SCENE 03------------------------------------------------------------
		float dist5;
		vec3 newPos5 = transform(transform(transform(pos + vec3(3.0, -1.0, 3.0), u_rotateY2), u_rotateX2), u_rotateZ2);	
		float bb5 = SDF_Sphere(newPos5, 1.1);
		if(bb5 < .015)
		{	
			float power = 10.0;
			dist5 = SDF_Mandelbulb(newPos5, power);
		}
		else
		{
			dist5 = bb5;
		}

		return dist5;
	}
}

vec3 backgroundColor() {
	int sn = sceneNum();
	float darken; 
	if(sn == 1)
	{
		darken = abs(cos(sin(8.0*f_uv.x*sin(u_time/12.0) + 8.0) + f_uv.y*2.0));
		darken *= abs(sin(cos(4.0*f_uv.y*2.0*sin(u_time/12.0) + 3.0) + f_uv.x*5.0));
	}
	else if(sn == 2)
	{
		darken = cos(48.0*length(f_uv - vec2(0.5, 0.5)) + sin(80.0*f_uv.x*-f_uv.y) + cos(50.0*-f_uv.x*f_uv.y) + sin(u_time));
	}
	else
	{
		darken = cos(length(f_uv - vec2(0.5, 0.5)));
		darken += (0.5 - length(vec2(0.25*f_uv.x + 0.25, f_uv.y) - vec2(0.5, 0.0)))/0.5;
	}
	
	darken = clamp(darken, 0.2, 1.0);

	vec3 a = vec3(0.5, 0.5, 0.5);
	vec3 b = vec3(0.5, 0.5, 0.5);
	vec3 c = vec3(2.0, 1.0, 1.0);
	vec3 d = vec3(0.5, 0.2, 0.25);
	float t = abs(sin(u_time/12.0));
	vec3 color = a + b*cos(6.28*(c*t + d));

	return darken*color;
}

// Compute the normal of an implicit surface using the gradient method
// Note: this method is slightly less accurate than sampling the scene at pos - epsilon, but because that is an expensive operation for the mandelbulb, we avoid it here
vec3 computeNormalFast( vec3 pos ) {
	vec2 point = vec2(0.0001, 0.0);
	float sampleAtPos = sceneMap(pos);
	vec3 normal = normalize(
			   vec3(sceneMap(pos + point.xyy) - sampleAtPos,
					sceneMap(pos + point.yxy) - sampleAtPos,
					sceneMap(pos + point.yyx) - sampleAtPos));
	return normal;
}

// Taken from a presentation by IQ: http://www.iquilezles.org/www/material/nvscene2008/rwwtt.pdf
float ComputeAO( vec3 pos, vec3 normal ) {
	float tStep = 0.0025;
	float diff = 0.0;
	float k = 34.0;
	for(float i = 1.0; i <= 5.0; i += 1.0) {
		vec3 sample = pos + (i * tStep) * normal;
		float dist = sceneMap( sample );
		diff += pow(0.5, i) * ((i * tStep) - dist);
	}
	return 1.0 - clamp(k * diff, 0.0, 0.9);
}

// Convert a fragment coordinate to normalized device coordinates
vec2 FragCoordToNDC( vec4 fragCoord ) {
	return 2.0 * vec2(fragCoord.x / u_resolution.x,
					  fragCoord.y / u_resolution.y) - 1.0;
}

// Return the direction of a ray cast through the given point in NDC
// This function includes the implementation of the camera; changes to the camera should be made here
vec3 Raycast( vec2 p_ndc, vec3 cameraPos ) {
	
	float len = 10.0;
	
	// Compute camera's frame of reference
	vec3 look = normalize(-cameraPos); // Assume we are looking towards (0, 0, 0)
	vec3 right = normalize(cross(look, vec3(0.0, 1.0, 0.0))); // 0, 1, 0 is the world up vector
	vec3 up = normalize(cross(right, look));
	
	float tanAlpha = tan(u_fovy / 2.0);
	vec3 V = up * len * tanAlpha;
	vec3 H = right * len * u_aspect * tanAlpha;
	
	// Convert x/y components of gl_FragCoord to NDC, then to a world space point
	vec3 point_World = p_ndc.x * H + p_ndc.y * V;
	
	// Return the direction
	return normalize(point_World - cameraPos);
}

// Referencing the following report: http://celarek.at/wp/wp-content/uploads/2014/05/realTimeFractalsReport.pdf
vec3 raymarchScene( vec3 origin, vec3 direction ) {
	// Basic Raymarching function: 
	float t = 0.01;
	vec3 result = vec3(0.0);
	
	for(int i = 1; i <= 750; i++) {
		float dist = sceneMap(origin + t * direction);
		
		if(abs(dist) < 0.0002) {
			result = vec3(t, float (i) / 1000.0, 1.0);
			break;
		} else if(t > T_MAX) {
			result = vec3(T_MAX, float (i) / 1000.0, 0.0);
			break;
		}
		
		#ifdef SPHERE_TRACING
			t += dist;
		#else
			t += 0.01;
		#endif
	}
	return result;
}

void main() {
	vec3 cameraPos = vec3(-3.5, 0, -3.5);
	vec2 point_NDC = FragCoordToNDC(gl_FragCoord);
	vec3 direction = Raycast(point_NDC, cameraPos);
	
	vec3 isect = raymarchScene( cameraPos, direction );
	vec3 isectPos = cameraPos + isect.x * direction;
	
	if(isect.z > 0.0) { // we did intersect with something
		vec3 normal = computeNormalFast( isectPos );
		normal = -normal;
		
		// Animated color using orbit traps: see the Mandelbulb SDF by IQ above
		vec3 trapColor;
		
		int sn = sceneNum();
		if(sn == 1)
		{
			trapColor = vec3(
				resColor.x-abs(sin((u_time) / 2.0 + isectPos.z) * 0.8), 
				resColor.y-(cos((u_time) / 8.0 + isectPos.y) * 0.5), 
				resColor.z+(cos((u_time) / 2.0 - isectPos.x) * 0.2));
		}
		else if(sn == 2)
		{
			trapColor = vec3(resColor) + vec3(0.4);
		}
		else
		{
			trapColor = vec3(0.0, resColor.y + 0.4, resColor.z + 0.4);
		}
		
		trapColor.r *= 0.01;
		
		// Two methods to compute AO:
		// float ao = ComputeAO(isectPos, normal);
		float fakeAO = 1.0 - clamp(isect.y, 0.0, 0.9);
		
		// gl_FragColor = ao * vec4(trapColor, 1);
		gl_FragColor = fakeAO * vec4(trapColor, 1);
	} else {
		gl_FragColor = vec4(backgroundColor(), 1);
	}
}
