
<script type="module">

	import Entry from "../Entry.svelte";
  import { onMount } from 'svelte';
  
  let p, q, r;
  // let setVars;

	import * as THREE from 'three';
	import { OrbitControls } from '$lib/OrbitControls.js';
	import { ImprovedNoise } from '$lib/ImprovedNoise.js';
	import { Stats } from '$lib/stats.module.js';

	let container, stats;

	let camera, controls, scene, renderer;

	let mesh, texture;

	let res = 1000;

	const worldWidth = 256, worldDepth = 256,
		worldHalfWidth = worldWidth / 2, worldHalfDepth = worldDepth / 2;

	let helper;

	const raycaster = new THREE.Raycaster();
	const pointer = new THREE.Vector2();

	onMount(() => {
		function setVars() {
      p = document.getElementById("p-entry").value;
      q = document.getElementById("q-entry").value;
      r = document.getElementById("r-entry").value;
    };

    setVars();
    document.addEventListener("input", setVars);

		init();
		animate();

		function init() {

			container = document.getElementById( 'container' );
			container.innerHTML = '';

			renderer = new THREE.WebGLRenderer( { antialias: true } );
			renderer.setPixelRatio( window.devicePixelRatio );
			renderer.setSize( window.innerWidth, window.innerHeight );
			container.appendChild( renderer.domElement );

			scene = new THREE.Scene();
			scene.background = new THREE.Color( 0xbfd1e5 );

			camera = new THREE.PerspectiveCamera( 60, window.innerWidth / window.innerHeight, 10, 20000 );

			controls = new OrbitControls( camera, renderer.domElement );
			controls.minDistance = 1000;
			controls.maxDistance = 10000;
			// controls.maxPolarAngle = Math.PI / 2;

			//

			const data = generateHeight( worldWidth, worldDepth );

			controls.target.y = data[ worldHalfWidth + worldHalfDepth * worldWidth ] + 500;
			camera.position.y = controls.target.y + 100;
			camera.position.x = 5000;
			camera.position.z = 500;
			controls.update();

			const geometry = new THREE.PlaneGeometry( res, res, worldWidth - 1, worldDepth - 1 );
			geometry.rotateX( - Math.PI / 2 );
			geometry.translate(res/2, 0, res/2)

			const vertices = geometry.attributes.position.array;

			for ( let i = 0, j = 0, l = vertices.length; i < l; i ++, j += 3 ) {

				vertices[ j + 1 ] = data[ i ] * 10;

			}

			//

			texture = new THREE.CanvasTexture( generateTexture( data, worldWidth, worldDepth ) );
			texture.wrapS = THREE.ClampToEdgeWrapping;
			texture.wrapT = THREE.ClampToEdgeWrapping;
			texture.colorSpace = THREE.SRGBColorSpace;

			let planeMat = new THREE.MeshBasicMaterial({
				color: 0xaa00aa, 
				// side: THREE.SingleSide, wireframe: true} 
				wireframe: true
			});

			mesh = new THREE.Mesh( geometry, new THREE.MeshBasicMaterial( { map: texture } ) );
			scene.add( mesh );

			// Add x-y plane
			const xy_geometry = new THREE.PlaneGeometry( res, res, worldWidth - 1, worldDepth - 1 );
			xy_geometry.translate(res/2, res/2, 0)

			const xy_vertices = geometry.attributes.position.array;
			let xy_plane = new THREE.Mesh( xy_geometry, planeMat );
			scene.add( xy_plane );

			// Add x-z plane
			const xz_geometry = new THREE.PlaneGeometry( res, res, worldWidth - 1, worldDepth - 1 );
			xz_geometry.rotateX( - Math.PI / 2 );
			xz_geometry.translate(res/2, 0, res/2)

			const xz_vertices = geometry.attributes.position.array;
			let xz_plane = new THREE.Mesh( xz_geometry, planeMat );
			scene.add( xz_plane );

			// Add y-z plane
			const yz_geometry = new THREE.PlaneGeometry( res, res, worldWidth - 1, worldDepth - 1 );
			yz_geometry.rotateY( + Math.PI / 2 );
			yz_geometry.translate(0, res/2, res/2)

			const yz_vertices = geometry.attributes.position.array;
			let yz_plane = new THREE.Mesh( yz_geometry, planeMat );
			scene.add( yz_plane );


			// const geometryHelper = new THREE.ConeGeometry( 20, 100, 3 );
			// geometryHelper.translate( 0, 50, 0 );
			// geometryHelper.rotateX( Math.PI / 2 );
			// helper = new THREE.Mesh( geometryHelper, new THREE.MeshNormalMaterial() );
			// scene.add( helper );

			// container.addEventListener( 'pointermove', onPointerMove );

			// stats = new Stats();
			// container.appendChild( stats.dom );

			//

			window.addEventListener( 'resize', onWindowResize );

		}

		function onWindowResize() {

			camera.aspect = window.innerWidth / window.innerHeight;
			camera.updateProjectionMatrix();

			renderer.setSize( window.innerWidth, window.innerHeight );

		}

		function generateHeight( width, height ) {

			const size = width * height, data = new Uint8Array( size ),
				perlin = new ImprovedNoise(), z = Math.random() * 100;

			let quality = 1;

			//for ( let j = 0; j < 4; j ++ ) {

				//for ( let i = 0; i < size; i ++ ) {

					//const x = i % width, y = ~ ~ ( i / width );
					//data[ i ] += Math.abs( perlin.noise( x / quality, y / quality, z ) * quality * 1.75 );

				//}

				//quality *= 5;

			//}

			console.log(data);
			return data;

		}

		function generateTexture( data, width, height ) {

			// bake lighting into texture

			let context, image, imageData, shade;

			const vector3 = new THREE.Vector3( 0, 0, 0 );

			const sun = new THREE.Vector3( 1, 1, 1 );
			sun.normalize();

			const canvas = document.createElement( 'canvas' );
			canvas.width = width;
			canvas.height = height;

			context = canvas.getContext( '2d' );
			context.fillStyle = '#000';
			context.fillRect( 0, 0, width, height );

			image = context.getImageData( 0, 0, canvas.width, canvas.height );
			imageData = image.data;

			for ( let i = 0, j = 0, l = imageData.length; i < l; i += 4, j ++ ) {

				vector3.x = data[ j - 2 ] - data[ j + 2 ];
				vector3.y = 2;
				vector3.z = data[ j - width * 2 ] - data[ j + width * 2 ];
				vector3.normalize();

				shade = vector3.dot( sun );

				imageData[ i ] = ( 96 + shade * 128 ) * ( 0.5 + data[ j ] * 0.007 );
				imageData[ i + 1 ] = ( 32 + shade * 96 ) * ( 0.5 + data[ j ] * 0.007 );
				imageData[ i + 2 ] = ( shade * 96 ) * ( 0.5 + data[ j ] * 0.007 );

			}

			context.putImageData( image, 0, 0 );

			// Scaled 4x

			const canvasScaled = document.createElement( 'canvas' );
			canvasScaled.width = width * 4;
			canvasScaled.height = height * 4;

			context = canvasScaled.getContext( '2d' );
			context.scale( 4, 4 );
			context.drawImage( canvas, 0, 0 );

			image = context.getImageData( 0, 0, canvasScaled.width, canvasScaled.height );
			imageData = image.data;

			for ( let i = 0, l = imageData.length; i < l; i += 4 ) {

				const v = ~ ~ ( Math.random() * 5 );

				imageData[ i ] += v;
				imageData[ i + 1 ] += v;
				imageData[ i + 2 ] += v;

			}

			context.putImageData( image, 0, 0 );

			return canvasScaled;

		}

		//

		function animate() {

			requestAnimationFrame( animate );

			render();
			// stats.update();

		}

		function render() {

			renderer.render( scene, camera );

		}

		function onPointerMove( event ) {

			pointer.x = ( event.clientX / renderer.domElement.clientWidth ) * 2 - 1;
			pointer.y = - ( event.clientY / renderer.domElement.clientHeight ) * 2 + 1;
			raycaster.setFromCamera( pointer, camera );

			// See if the ray from the camera into the world hits one of our meshes
			const intersects = raycaster.intersectObject( mesh );

			// Toggle rotation bool for meshes that we clicked
			if ( intersects.length > 0 ) {

				helper.position.set( 0, 0, 0 );
				helper.lookAt( intersects[ 0 ].face.normal );

				helper.position.copy( intersects[ 0 ].point );

			}

		}
	})
</script>


<h1>Mean Minimal Spanning Tree Calculator</h1>

<row style="display: flex; flex-flow: row; justify-content: space-around;">
<column>
<Entry name="p"/>
<Entry name="q"/>
<Entry name="r"/>
</column>
<div>V = </div>
<column>
0 {p} {q} {r}<br>{p} 0 {r} {q}<br>{q} {r} 0 {p}<br>{r} {q} {p} 0
</column>
</row>


<div id="container"></div>
<div id="info"><a href="https://threejs.org" target="_blank" rel="noopener">three.js</a> - webgl terrain raycasting demo</div>
