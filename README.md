# Parallelization of Physical Particle System Simulation
## About
本次實作主要關注於如何將平行化技術應用在動畫特效上，期望能夠優化動畫粒子的運算，為觀眾帶來更好的視覺饗宴，同時我們會針對此一問題，利用C++與OpenGL設計一套物理粒子模擬系統，分別為sequential與parallel的版本，為了讓優化效果更加顯著，我們藉由提升物理動畫品質來提升運算量，讓sequential版本的程式到達效能瓶頸。而後我們會運用課堂上所介紹的平行化技術，將其套用在動畫繪製上，而動畫品質的好壞，則預計採用FPS作為指標。實驗中我們會觀察並記錄動畫的FPS變化，比較彼此之間所帶來的效能優化幅度，並從中討論歸納各種環境因素對其造成的影響，最後做出結論統整。
## Related Works
- [OpenGL 互動式遊戲](https://github.com/MarkMoHR/OpenglGame)
- [OpenGL搭配GPU加速粒子系統模擬](https://github.com/MauriceGit/Partikel_accelleration_on_GPU)
- [OpenGL粒子模擬器](https://github.com/Syntaf/ParticleSimulator)
- [OpenGL空間碰撞檢測](https://github.com/ttvd/spatial-collision-datastructures)
- [2D粒子碰撞模擬](https://github.com/nicolasruan/elastic-collisions-2d)
## Language and Library
- C++
- OpenGL
- Pthread
- OpenMP
- CUDA
## Job Distribution
||Brute Force|Sort & Sweep|Uniform Grid|Hierarchical Grid|Octree|Loose Octree|KDTree|
|-|-|-|-|-|-|-|-|
|In Charge|JoeHuang99|JoeHuang99|JoeHuang99|niangao19|eesss34690|eesss34690|niangao19|
