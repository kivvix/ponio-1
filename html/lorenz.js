import * as THREE from 'three';

const camera_radius = 100.;
const scene = new THREE.Scene();
const camera = new THREE.PerspectiveCamera( 45, window.innerWidth / window.innerHeight, 1, 500 );
camera.position.set( 0, 0, camera_radius );
camera.lookAt( 0, 0, 0 );

const renderer = new THREE.WebGLRenderer({antialias: true, alpha: true});
renderer.setSize( window.innerWidth, window.innerHeight );
renderer.setAnimationLoop( animate );
let canvas = renderer.domElement;
canvas.id = "lorenz-canvas";
document.body.appendChild( canvas );

let tn = 0;

function lorenz(y0, dy)
{
  let sigma = 10.;
  let rho = 20.;
  let beta = 8./3.;

  dy.x = sigma*(y0.y - y0.x);
  dy.y = rho*y0.x - y0.y - y0.x*y0.z;
  dy.z = y0.x*y0.y - beta*y0.z;
}

class runge_kutta
{
  dt = 0.005;
  y0 = new THREE.Vector3( 1., 1., 1. );
  y1 = new THREE.Vector3( 1., 1., 1. );
  points = [];

  constructor(color, n_stages) {
    this.material = new THREE.LineBasicMaterial( { color: color } );
    this.stages = Array(n_stages).fill(new THREE.Vector3( 0., 0., 0. ));
  }

  next() {
    this.points.push(new THREE.Vector3(this.y0.x, this.y0.y, this.y0.z));
    this.increment();

    this.y0 = this.y1;
  }

  geometry() {
    if (this.points.length < 2 )
    {
      return new THREE.BufferGeometry().setFromPoints( this.points );
    }
    return new THREE.BufferGeometry().setFromPoints([this.points[this.points.length-2], this.points[this.points.length-1]] );
  }

  add_tO_scene(scene) {
    const geometry = this.geometry();
    const line = new THREE.Line( geometry, this.material );

    scene.add(line);
  }
};

let euler = new runge_kutta(0x079992, 1);
euler.increment = function () {
  lorenz(this.y0, this.stages[0]);

  this.y1.x = this.y0.x + this.dt*this.stages[0].x;
  this.y1.y = this.y0.y + this.dt*this.stages[0].y;
  this.y1.z = this.y0.z + this.dt*this.stages[0].z;
};

let rk3 = new runge_kutta(0xb71540, 3);
rk3.increment = function() {
  lorenz(this.y0, this.stages[0]);

  this.y1.x = this.y0.x + this.dt*this.stages[0].x;
  this.y1.y = this.y0.y + this.dt*this.stages[0].y;
  this.y1.z = this.y0.z + this.dt*this.stages[0].z;

  lorenz(this.y1, this.stages[1]);

  this.y1.x = this.y0.x + 0.25*this.dt*(this.stages[0].x + this.stages[1].x);
  this.y1.y = this.y0.y + 0.25*this.dt*(this.stages[0].y + this.stages[1].y);
  this.y1.z = this.y0.z + 0.25*this.dt*(this.stages[0].z + this.stages[1].z);

  lorenz(this.y1, this.stages[2]);

  this.y1.x = this.y0.x + this.dt*(this.stages[0].x/6. + this.stages[1].x/6. + 2./3.*this.stages[2].x);
  this.y1.y = this.y0.y + this.dt*(this.stages[0].y/6. + this.stages[1].y/6. + 2./3.*this.stages[2].y);
  this.y1.z = this.y0.z + this.dt*(this.stages[0].z/6. + this.stages[1].z/6. + 2./3.*this.stages[2].z);
};

function animate() {

  euler.next();
  rk3.next();
  euler.add_tO_scene(scene);
  rk3.add_tO_scene(scene);

  tn += 0.005;

  camera.position.z = camera_radius*Math.sin(tn);
  camera.position.x = camera_radius*Math.cos(tn);

  camera.lookAt( 0, 0, 0 );


  renderer.render( scene, camera );
}
