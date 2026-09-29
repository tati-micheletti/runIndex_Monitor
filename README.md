After knitting and rendering `runIndex_Monitor.Rmd`, overwrite this file with a symlink to `runIndex_Monitor.md`.

E.g., on Linux/macOS from terminal, in the directory that contains the module:

```bash
cd runIndex_Monitor
rm README.md && ln -s runIndex_Monitor.md README.md
```
