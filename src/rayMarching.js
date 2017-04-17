const THREE = require('three');
const EffectComposer = require('three-effectcomposer')(THREE)

import {PROXY_BUFFER_SIZE} from './proxy_geometry'

export default function RayMarcher(renderer, scene, camera) {
    var composer = new EffectComposer(renderer);
    var shaderPass = new EffectComposer.ShaderPass({
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
            }
        },
        vertexShader: require('./glsl/pass-vert.glsl'),
        fragmentShader: require('./glsl/rayMarch-frag.glsl')
        
    });
    shaderPass.renderToScreen = true;
    composer.addPass(shaderPass);

    return {
        render: function(buffer, clock) {
            shaderPass.uniforms["u_time"].value = clock.getElapsedTime();
            composer.render();
            // console.log(composer);
        }
    }
}