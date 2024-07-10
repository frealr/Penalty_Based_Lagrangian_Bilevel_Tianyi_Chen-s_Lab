# TFG-BLO. Algoritmo BLOCC

## Carpeta `hyperparam_opt` 
La carpeta `hyperparam_opt` contiene experimentos sobre la optimización de hiperparámetros de SVM. Se utiliza comúnmente para probar el rendimiento de los algoritmos BLO, tales como LV-HBA [Yao et al., 2023] y GA [Xu y Zhu, 2023]. Los detalles de sus experimentos también se pueden ver en:
- [LV-HBA](https://github.com/SUSTech-Optimization/LV-HBA)
- [Efficient Gradient Approximation Method](https://github.com/xsy786912649/Efficient-gradient-approximation-method-for-constrained-bilevel-optimization-problem)

### Contenido
- Tres algoritmos están en archivos py en la carpeta `algorithms`. Los algoritmos ours (BLOCC) y LV-HBA resuelven dos formulaciones diferentes. Cada una va asociada a un fichero, teniendo ours.py, ours_reformulacion.py, lv_hba.py y lv_hba_reformulacion.py
- Los conjuntos de datos incluyen el conjunto de datos de diabetes y el conjunto de datos fourclass.
- Por favor, ejecutar el notebook de Jupyter más actualizado. Los 4 últimos dígitos hacen referencia a la fecha de creación (dd/mm) de 2024.

## Carpeta `Transportation` 
La carpeta `Transportation` contiene experimentos sobre el problema de diseño de redes de transporte.

## `toy_example.ipynb`
El archivo `toy_example.ipynb` presenta 2 ejemplos sencillos para mostrar que el Algoritmo BLOCC funciona eficazmente.
