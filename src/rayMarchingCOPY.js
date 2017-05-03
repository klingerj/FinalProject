const THREE = require('three');
const EffectComposer = require('three-effectcomposer')(THREE)

export default function RayMarcher(renderer, scene, camera) {
	
	var target1 = new THREE.WebGLRenderTarget(window.innerWidth, window.innerHeight);
    var composer1 = new EffectComposer(renderer, target1);
    var shaderPass1 = new EffectComposer.ShaderPass({
        uniforms: {
            u_time: {
                type: 'f',
                value: 0
            },
            u_resolution: {
                type: 'v2',
                value: new THREE.Vector2(window.innerWidth, window.innerHeight)
            },
            u_fovy: {
                type: 'f',
                value: camera.fov
            },
            u_aspect: {
                type: 'f',
                value: camera.aspect
            },
            u_cwMat: {
                type: 'm4v',
                value: null
            },
            u_ccwMat: {
                type: 'm4v',
                value: null
            },
            u_northMat: {
                type: 'm4v',
                value: null
            },
            u_southMat: {
                type: 'm4v',
                value: null
            },
            u_westMat: {
                type: 'm4v',
                value: null
            },
            u_eastMat: {
                type: 'm4v',
                value: null
            },
            u_rotateX1: {
                type: 'm4v',
                value: null
            },
            u_rotateY1: {
                type: 'm4v',
                value: null
            },
            u_rotateZ1: {
                type: 'm4v',
                value: null
            },
            u_rotateX2: {
                type: 'm4v',
                value: null
            },
            u_rotateY2: {
                type: 'm4v',
                value: null
            },
            u_rotateZ2: {
                type: 'm4v',
                value: null
            }
        },
        vertexShader: require('./glsl/pass-vert.glsl'),
        fragmentShader: require('./glsl/firstPass-frag.glsl')
    });
    
    var composer2 = new EffectComposer(renderer);
    var shaderPass2 = new EffectComposer.ShaderPass({
        uniforms: {
            u_firstPass: {
            	type: 't',
            	value: null
            },
            u_previousFrame: {
            	type: 't',
            	value: null
            }
        },
        vertexShader: require('./glsl/pass-vert.glsl'),
        fragmentShader: require('./glsl/secondPass-frag.glsl')
    });
    shaderPass2.renderToScreen = true;
    shaderPass2.material.uniforms.u_firstPass.value = composer1.writeBuffer.texture;
    
    var target3 = new THREE.WebGLRenderTarget(window.innerWidth, window.innerHeight);
    var composer3 = new EffectComposer(renderer, target3);
    var shaderPass3 = new EffectComposer.ShaderPass({
        uniforms: {
            u_input: {
            	type: 't',
            	value: null
            }
        },
        vertexShader: require('./glsl/pass-vert.glsl'),
        fragmentShader: require('./glsl/thirdPass-frag.glsl')
    });
    shaderPass2.material.uniforms.u_previousFrame.value = composer3.writeBuffer.texture;
    shaderPass3.material.uniforms.u_input.value = composer1.writeBuffer.texture;
    
    composer1.addPass(shaderPass1);
    composer2.addPass(shaderPass2);
    composer3.addPass(shaderPass3);

    return {
        render: function(buffer, clock) {
            shaderPass1.uniforms["u_time"].value = clock.getElapsedTime();
            
            // Mandelbulb transformation uniforms
            var angle = clock.getElapsedTime() / (4.0 * 3.1415); 

            var cwMat = new THREE.Matrix4();
            cwMat.makeRotationY(angle);

            var ccwMat = new THREE.Matrix4();
            ccwMat.makeRotationY(-angle);

            var northMat = new THREE.Matrix4();
            northMat.makeRotationX(angle);

            var southMat = new THREE.Matrix4();
            southMat.makeRotationX(-angle);

            var eastMat = new THREE.Matrix4();
            eastMat.makeRotationZ(angle);

            var westMat = new THREE.Matrix4();
            westMat.makeRotationZ(-angle);

            var rotateX1 = new THREE.Matrix4();
            rotateX1.makeRotationX(3.1415 / 2.0);

            var rotateY1 = new THREE.Matrix4();
            rotateY1.makeRotationY(45.0 * 3.1415 / 180.0);

            var rotateZ1 = new THREE.Matrix4();
            rotateZ1.makeRotationZ((2.0 * clock.getElapsedTime()) * 3.1415 / 180.0);

            var rotateX2 = new THREE.Matrix4();
            rotateX2.makeRotationX(3.1415 / 2.0);

            var rotateY2 = new THREE.Matrix4();
            rotateY2.makeRotationY(-45.0 * 3.1415 / 180.0);

            var rotateZ2 = new THREE.Matrix4();
            rotateZ2.makeRotationY((2.0 * clock.getElapsedTime()) * 3.1415 / 180.0);
            
            shaderPass1.uniforms["u_cwMat"].value = cwMat;
            shaderPass1.uniforms["u_ccwMat"].value = ccwMat;
            shaderPass1.uniforms["u_northMat"].value = northMat;
            shaderPass1.uniforms["u_southMat"].value = southMat;
            shaderPass1.uniforms["u_westMat"].value = westMat;
            shaderPass1.uniforms["u_eastMat"].value = eastMat;
            shaderPass1.uniforms["u_rotateX1"].value = rotateX1;
            shaderPass1.uniforms["u_rotateY1"].value = rotateY1;
            shaderPass1.uniforms["u_rotateZ1"].value = rotateZ1;
            shaderPass1.uniforms["u_rotateX2"].value = rotateX2;
            shaderPass1.uniforms["u_rotateY2"].value = rotateY2;
            shaderPass1.uniforms["u_rotateZ2"].value = rotateZ2;
            
            composer1.render();
            composer2.render();
            composer3.render();
        }
    }
}