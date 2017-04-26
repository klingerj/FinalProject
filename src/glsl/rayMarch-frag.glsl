
#define MAX_GEOMETRY_COUNT 100
#define SPHERE_TRACING true
#define T_MAX 20.0

/* This is how I'm packing the data
struct geometry_t {
    vec3 position;
    float type;
};
*/
// uniform vec4 u_buffer[MAX_GEOMETRY_COUNT];
// uniform int u_count;

varying vec2 f_uv;

uniform float u_time;
uniform vec2 u_resolution;
uniform float u_fovy;
uniform float u_aspect;

vec4 resColor;

/***** Geometry SDF Functions
http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
							  *****/

float SDF_Sphere( vec3 pos, float radius ) {
	return length(pos) - radius;
}

//diagonal is the vector from the center of the box to the first quadrant corner
float boxSDF(vec3 point, vec3 diagonal) {
	vec3 d = abs(point) - diagonal;
  	return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

float SDF_Mandlebulb( vec3 p , float manPower)
{
	vec3 w = p;
    float m = dot(w,w);

    vec4 trap = vec4(abs(w),m);
    float dz = 1.0;
    
    
    for( int i=0; i<4; i++ )
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

    return 0.25*log(m)*sqrt(m)/dz;
}

//Operators:

float intersection(float d1, float d2)
{
    return max(d1,d2);
}

float subtraction( float d1, float d2 )
{
    return max(-d1,d2);
}

float un(float d1, float d2)
{
    return min(d1,d2);
}

//returns transformed point based on rotation and translation matrix of shape
vec3 transform(vec3 point, mat4 trans)
{
	//columns of the rotation matrix transpose
	vec3 col1 = vec3(trans[0][0], trans[1][0], trans[2][0]);
	vec3 col2 = vec3(trans[0][1], trans[1][1], trans[2][1]);
	vec3 col3 = vec3(trans[0][2], trans[1][2], trans[2][2]);

	mat3 rotTranspose = mat3(col1, col2, col3);

	vec3 col4 = -1.0*rotTranspose*vec3(trans[3]);

	mat4 newTrans = mat4(vec4(col1, 0.0), vec4(col2, 0.0), vec4(col3, 0.0), vec4(col4, 1.0));

	return vec3(newTrans * vec4(point, 1.0));
}

// Return the distance of the closest object in the scene
float sceneMap( vec3 pos ) {
	return SDF_Sphere( pos, 1.0 );
}

float sceneMap2( vec3 pos ){

	float angle = u_time/(2.0*3.1415);
	mat4 cwMat = mat4(1.0); //transform for moving clockwise
	cwMat[0][0] = cos(angle); cwMat[0][2] = -sin(angle); cwMat[2][0] = sin(angle); cwMat[2][2] = cos(angle); //rotating about y-axis, based on utime
	mat4 ccwMat = mat4(1.0); //transform for moving counterclockwise
	ccwMat[0][0] = cos(-angle); ccwMat[0][2] = -sin(-angle); ccwMat[2][0] = sin(-angle); ccwMat[2][2] = cos(-angle); //rotating about y-axis, based on utime

	mat4 northMat = mat4(1.0); 
	northMat[1][1] = cos(angle); northMat[1][2] = sin(angle); northMat[2][1] = -sin(angle); northMat[2][2] = cos(angle); //rotating about x-axis, based on utime
	mat4 southMat = mat4(1.0); 
	southMat[1][1] = cos(-angle); southMat[1][2] = sin(-angle); southMat[2][1] = -sin(-angle); southMat[2][2] = cos(-angle); //rotating about x-axis, based on utime
	mat4 westMat = mat4(1.0); 
	westMat[0][0] = cos(-angle); westMat[0][1] = sin(-angle); westMat[1][0] = -sin(-angle); westMat[1][1] = cos(-angle); //rotating about z-axis, based on utime
	mat4 eastMat = mat4(1.0); 
	eastMat[0][0] = cos(angle); eastMat[0][1] = sin(angle); eastMat[1][0] = -sin(angle); eastMat[1][1] = cos(angle); //rotating about z-axis, based on utime

	// vec3 newPos1 = transform(pos + vec3(0, 1.5, 0), cwMat);
	// vec3 newPos2 = transform(transform(pos + vec3(1.5, 0, 0), ccwMat), eastMat);
	// vec3 newPos3 = transform(transform(pos + vec3(-1.5, 0, 0), ccwMat), westMat);
	// vec3 newPos4 = transform(transform(pos + vec3(0, 0, 1.5), ccwMat), northMat);
	// vec3 newPos5 = transform(transform(pos + vec3(0, 0, -1.5), ccwMat), southMat);
	// vec3 newPos6 = transform(pos + vec3(2, -1.5, 2), cwMat);
	// vec3 newPos7 = transform(pos + vec3(-2, -1.5, -2), cwMat);
	// vec3 newPos8 = transform(pos + vec3(-2, -1.5, 2), cwMat);
	// vec3 newPos9 = transform(pos + vec3(2, -1.5, -2), cwMat);
	
	// float dist1;
	// float bb1 = boxSDF(newPos1, vec3(1.5,1.5,1.5));
	// if(bb1 < .015)
	// {
	// 	float power = 10.0;//12.0 + abs(sin(u_time/4.0))*40.0;
	// 	dist1 = SDF_Mandlebulb(newPos1, power);
	// }
	// else
	// {
	// 	dist1 = bb1;
	// }

	// float dist2;
	// float bb2 = boxSDF(newPos2, vec3(1.1,1.1,1.1));
	// if(bb2 < .015)
	// {
	// 	dist2 = SDF_Mandlebulb(newPos2, 16.0);
	// }
	// else
	// {
	// 	dist2 = bb2;
	// }

	// float dist3;
	// float bb3 = boxSDF(newPos3, vec3(1.1,1.1,1.1));
	// if(bb3 < .015)
	// {
	// 	dist3 = SDF_Mandlebulb(newPos3, 16.0);
	// }
	// else
	// {
	// 	dist3 = bb3;
	// }

	// float dist4;
	// float bb4 = boxSDF(newPos4, vec3(1.1,1.1,1.1));
	// if(bb4 < .015)
	// {
	// 	dist4 = SDF_Mandlebulb(newPos4, 16.0);
	// }
	// else
	// {
	// 	dist4 = bb4;
	// }

	// float dist5;
	// float bb5 = boxSDF(newPos5, vec3(1.1,1.1,1.1));
	// if(bb5 < .015)
	// {
	// 	dist5 = SDF_Mandlebulb(newPos5, 16.0);
	// }
	// else
	// {
	// 	dist5 = bb5;
	// }

	//float man2 = SDF_Mandlebulb(newPos2, 16.0);
	//float man3 = SDF_Mandlebulb(newPos3, 16.0);
	//float man4 = SDF_Mandlebulb(newPos4, 16.0);
	//float man5 = SDF_Mandlebulb(newPos5, 16.0);
	//float man6 = SDF_Mandlebulb(newPos6, 24.0);
	//float man7 = SDF_Mandlebulb(newPos7, 24.0);
	//float man8 = SDF_Mandlebulb(newPos8, 24.0);
	//float man9 = SDF_Mandlebulb(newPos9, 24.0);
	
	//return un(dist1, un(dist2, un(dist3, un(dist4, dist5))));
	//return un(man1, un(man2, un(man3, un(man4, un(man5, un(man6, un(man7, un(man8, man9))))))));

	float dist1;
	vec3 newPos1 = transform(transform(pos + vec3(cos((u_time+4.0)/8.0)*4.0, 0, sin(u_time/7.0)*3.5), cwMat), northMat);	
	float bb1 = boxSDF(newPos1, vec3(1.1,1.1,1.1));
	if(bb1 < .015)
	{
		float power = 12.0;
		dist1 = SDF_Mandlebulb(newPos1, power);
	}
	else
	{
		dist1 = bb1;
	}

	float dist2;
	vec3 newPos2 = transform(transform(pos + vec3(cos((u_time+50.0)/10.0)*2.0, 1, sin((u_time+30.0)/6.0)*2.5), ccwMat), eastMat);
	float bb2 = boxSDF(newPos2, vec3(1.1,1.1,1.1));
	if(bb2 < .015)
	{
		
		float power = 12.0;
		dist2 = SDF_Mandlebulb(newPos2, power);
	}
	else
	{
		dist2 = bb2;
	}

	float dist3;
	vec3 newPos3 = transform(transform(pos + vec3(sin((u_time)/16.0)*3.0, -1, cos((u_time+75.0)/3.0)*2.0), cwMat), westMat);
	float bb3 = boxSDF(newPos3, vec3(1.1,1.1,1.1));
	if(bb3 < .015)
	{
		
		float power = 12.0;
		dist3 = SDF_Mandlebulb(newPos3, power);
	}
	else
	{
		dist3 = bb3;
	}

	return un(dist1, un(dist2, dist3));

	// float dist4;
	// float bb4 = boxSDF(newPos4, vec3(1.2,1.2,1.2));
	// if(bb4 < .015)
	// {
	// 	dist4 = SDF_Mandlebulb(newPos4, 16.0);
	// }
	// else
	// {
	// 	dist4 = bb4;
	// }

	// float dist5;
	// float bb5 = boxSDF(newPos5, vec3(1.2,1.2,1.2));
	// if(bb5 < .015)
	// {
	// 	dist5 = SDF_Mandlebulb(newPos5, 16.0);
	// }
	// else
	// {
	// 	dist5 = bb5;
	// }
	
	//float man2 = SDF_Mandlebulb(newPos2, 16.0);
	//float man3 = SDF_Mandlebulb(newPos3, 16.0);
	//float man4 = SDF_Mandlebulb(newPos4, 16.0);
	//float man5 = SDF_Mandlebulb(newPos5, 16.0);
	//float man6 = SDF_Mandlebulb(newPos6, 24.0);
	//float man7 = SDF_Mandlebulb(newPos7, 24.0);
	//float man8 = SDF_Mandlebulb(newPos8, 24.0);
	//float man9 = SDF_Mandlebulb(newPos9, 24.0);
	// return un(dist1, un(dist2, un(dist3, un(dist4, dist5))));
	
	// dist1 = SDF_Mandlebulb(newPos1, 12.0);
	//return dist1;
	//return un(man1, un(man2, un(man3, un(man4, un(man5, un(man6, un(man7, un(man8, man9))))))));
}

// Compute the normal of an implicit surface using the gradient method
vec3 computeNormal( vec3 pos ) {
	vec2 point = vec2(0.0001, 0.0);
	vec3 normal = normalize(
			   vec3(sceneMap2(pos + point.xyy) - sceneMap2(pos - point.xyy),
					sceneMap2(pos + point.yxy) - sceneMap2(pos - point.yxy),
					sceneMap2(pos + point.yyx) - sceneMap2(pos - point.yyx)));
	return normal;
}

// Check for intersection with the scene for increasing t-values
vec2 raymarchScene( vec3 origin, vec3 direction ) {
	float dist;
	float t = 0.01;
	for(int i = 0; i < 500; i++) {
		float dist = sceneMap2(origin + t * direction);
		if(dist < 0.0001) {
			return vec2(t, 1.0); // intersection
		} else if(t > T_MAX) {
			break;
		}
		#ifdef SPHERE_TRACING
			t += dist;
		#else
			t += 0.01;
		#endif
	}
	return vec2(0.0, -1.0); // no intersection
}


float SpecHighlight( vec3 toCam, vec3 toLight, vec3 normal) {
	float dot = dot(normalize(toCam + toLight), normal);
	return max(dot * dot * dot * dot * dot * dot * dot * dot, 0.0);
}

// Presentation by IQ: http://www.iquilezles.org/www/material/nvscene2008/rwwtt.pdf
float ComputeAO( vec3 pos, vec3 normal ) {
	float tStep = 0.0025;
	float t = 0.0;
	float ao = 1.0;
	float diff = 0.0;
	float k = 72.0;
	for(int i = 0; i < 5; i++) {
		vec3 sample = pos + t * normal;
		float dist = sceneMap2( sample );
		diff += pow(0.5, float (i)) * (t - dist);
		t += tStep;
	}
	ao -= clamp(k * diff, 0.0, 1.0);
	return ao;
}



void main() {
	
	/** Raycasting **/
	
	// Convering gl_FragCoord to normalized device coordinates: http://www.txutxi.com/?p=182
	vec2 point_NDC = 2.0 * vec2(gl_FragCoord.x / u_resolution.x,
								gl_FragCoord.y / u_resolution.y) - 1.0;

	vec3 cameraPos = vec3(-3.5, 0, -3.5);
	//vec3 cameraPos = vec3(1, -4, 2);
	//vec3 cameraPos = vec3(1, 0, 1);
	
	// Circle the origin (0, 0, 0)
	// cameraPos.x = sin(u_time) * 10.0;
	// cameraPos.z = cos(u_time) * 10.0;
	
	float len = 10.0; // assume the reference point is at 0, 0, 0
	
	
	// Compute camera's frame of reference
	vec3 look = normalize(-cameraPos);
	vec3 right = normalize(cross(look, vec3(0.0, 1.0, 0.0))); // 0, 1, 0 is the world up vector
	vec3 up = normalize(cross(right, look));
	
	float tanAlpha = tan(u_fovy / 2.0);
	vec3 V = up * len * tanAlpha;
	vec3 H = right * len * u_aspect * tanAlpha;
	
	// Convert x/y components of gl_FragCoord to NDC, then to a world space point
	vec3 point_World = point_NDC.x * H + point_NDC.y * V;
	
	// Perform the raymarch
	vec3 direction = normalize(point_World - cameraPos);
	vec2 isect = raymarchScene( cameraPos, direction );
	vec3 isectPos = cameraPos + isect.x * direction;
	
	/** Shading and lighting **/
	
	if(isect.y > 0.0) { // we did intersect with something
		vec3 normal = computeNormal( isectPos );
		
		// Lighting
		vec3 baseMaterial = vec3(0.2, 0.2, 0.2);
		vec3 trapColor = vec3(resColor.x+(sin((u_time+isectPos.x)/2.0)*0.2), resColor.y-(cos((u_time+isectPos.y)/8.0)*0.5), resColor.z+(cos((u_time-isectPos.z)/2.0)*0.8));
		vec3 sun = vec3(0.5, 0.4, 0.3) * 12.0;
		vec3 sunPos = vec3(5.0, 5.0, 0.0);
		
		vec3 toSun = normalize(sunPos - isectPos);
		
		normal = -normal;
		
		// Visibility test
		// vec2 shadowTest = raymarchScene( isectPos, toSun );
		// float vis = 1.0;
		
		// if(shadowTest.y > 0.0) { // something is blocking this point
		// 	vis = 0.0;
		// }
		
		
		
		// Phong-ish shading for now
		float spec = SpecHighlight( -look, toSun, normal);
		
		float sunDot = clamp(dot( normal, toSun ), 0.0, 1.0);
		float ao = ComputeAO(isectPos, normal);
		
		// Apply lambertian shading - for now
		gl_FragColor = /*vis * */ao * vec4(((1.0 - spec) * trapColor * baseMaterial * sun /** vec3(sunDot)*/ + spec * vec3(0.1)), 1);
		//gl_FragColor = vec4(clamp(normal.x, 0.1, 0.9), -normal.y, normal.z, 1);
	} else {
		// Background color
		gl_FragColor = vec4(0.5, 0.5, 0.5, 1);
	}
}
