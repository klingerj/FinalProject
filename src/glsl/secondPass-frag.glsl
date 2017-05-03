/* 
	Code written by Joseph Klinger and Tabatha Hickman
	University of Pennsylvania
	CIS 700 Instructor: Rachel Hwang
	May 2017
*/

#define SPHERE_TRACING true
#define T_MAX 12.0
#define MOTION_BLUR true

varying vec2 f_uv;

uniform sampler2D u_firstPass;
uniform sampler2D u_previousFrame;

// The sole purpose of this pass is to compute motion blur between the current and previous frames
void main() {
	#ifdef MOTION_BLUR
		for(float w = 0.0; w <= 1.0; w += 0.2) {
			gl_FragColor += (1.0 - w) * texture2D(u_firstPass, f_uv) + w * texture2D(u_previousFrame, f_uv);
		}
		gl_FragColor /= 1.0 / 0.2; // divide by the number of samples (1 / stepSize)
	#else
		gl_FragColor = texture2D(u_firstPass, f_uv);
	#endif
}
