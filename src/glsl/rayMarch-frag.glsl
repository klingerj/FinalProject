
#define MAX_GEOMETRY_COUNT 100
#define SPHERE_TRACING true
#define T_MAX 40.0

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

/***** Geometry SDF Functions
http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm
							  *****/

float SDF_Sphere( vec3 pos, float radius ) {
	return length(pos) - radius;
}

// Return the distance of the closest object in the scene
float sceneMap( vec3 pos ) {
	return SDF_Sphere( pos, 1.0 );
}

// Compute the normal of an implicit surface using the gradient method
vec3 computeNormal( vec3 pos ) {
	vec2 point = vec2(0.0001, 0.0);
	vec3 normal = normalize(
			   vec3(sceneMap(pos + point.xyy) - sceneMap(pos - point.xyy),
					sceneMap(pos + point.yxy) - sceneMap(pos - point.yxy),
					sceneMap(pos + point.yyx) - sceneMap(pos - point.yyx)));
	return normal;
}

// Check for intersection with the scene for increasing t-values
vec2 raymarchScene( vec3 origin, vec3 direction ) {
	float dist;
	float t = 0.01;
	for(int i = 0; i < 500; i++) {
		float dist = sceneMap(origin + t * direction);
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






void main() {
	
	/** Raycasting **/
	
	// Convering gl_FragCoord to normalized device coordinates: http://www.txutxi.com/?p=182
	vec2 point_NDC = 2.0 * vec2(gl_FragCoord.x / u_resolution.x,
								gl_FragCoord.y / u_resolution.y) - 1.0;

	vec3 cameraPos = vec3(0.0, 0.0, 0.0);
	
	// Circle the origin (0, 0, 0)
	cameraPos.x = sin(u_time) * 10.0;
	cameraPos.z = cos(u_time) * 10.0;
	
	float len = 20.0; // assume the reference point is at 0, 0, 0
	
	
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
		vec3 baseMaterial = vec3(0.2);
		vec3 sun = vec3(0.5, 0.4, 0.3) * 12.0;
		vec3 sunPos = vec3(5.0, 5.0, 0.0);
		float sunDot = clamp(dot( -normal, normalize(sunPos - isectPos) ), 0.0, 1.0);
		
		// Apply lambertian shading - for now
		gl_FragColor = vec4( baseMaterial * sun * vec3(sunDot), 1 );
	} else {
		// Background color
		gl_FragColor = vec4(f_uv, 0, 1);
	}
}
