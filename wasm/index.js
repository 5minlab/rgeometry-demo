// Note that a dynamic `import` statement here is required due to
// webpack/webpack#6615, but in theory `import { greet } from './pkg';`
// will work here one day as well!
const rust = import('./node/asm.js');

rust
  .then(m => {
      let sx = m.Simplical.new();

      for (let i = 0; i < 10; i++) {
          const sx0 = m.Simplical.from_rect(i+0.1, i+0.1, 2, 2, 0);
          const sx1 = sx.union(sx0);

          sx.free();
          sx0.free();
          sx = sx1;
      }

      const triangulated = m.Triangulated.from(sx);
      const visibility = triangulated.visibility(0, 0);
      console.log(visibility);
      triangulated.free();
  })
  .catch(console.error);
