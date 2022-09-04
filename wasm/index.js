// Note that a dynamic `import` statement here is required due to
// webpack/webpack#6615, but in theory `import { greet } from './pkg';`
// will work here one day as well!
const rust = import('./pkg/asm.js');

rust
  .then(m => {
      let sx = m.simplical_new();

      for (let i = 0; i < 10; i++) {
          const p = m.gen_points_cube([i+0.1, i+0.1], 2);
          const sx0 = m.simplical_from_polygon(p);
          const sx1 = m.simplical_union(sx, sx0);

          m.simplical_free(sx);
          m.simplical_free(sx0);
          sx = sx1;
      }

      const simplices = m.simplical_simplices(sx);
      console.log(simplices);

      const triangulated = m.Triangulated.from(sx);
      const visibility = triangulated.visibility(0, 0);
      console.log(visibility);
      triangulated.free();
  })
  .catch(console.error);
