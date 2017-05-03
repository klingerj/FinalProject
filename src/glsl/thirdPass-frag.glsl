/* 
	Code written by Joseph Klinger and Tabatha Hickman
	University of Pennsylvania
	CIS 700 Instructor: Rachel Hwang
	May 2017
*/

varying vec2 f_uv;

uniform sampler2D u_input;

void main() {
	gl_FragColor = texture2D(u_input, f_uv);
}
