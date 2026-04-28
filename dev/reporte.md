# Teoría Espectral de Gráficas

## Portada
> **Proyecto**: Análisis espectral de la estructura y dinámica de redes financieras.

**Autor**: Dylan Ramírez Hernández

## Índice
1. [Introducción](#1-introducción)
2. [Marco Teórico](#2-marco-teórico)

   2.1 [Mercado financiero como red compleja](#21-mercado-financiero-como-red-compleja)  
   2.2 [Operadores Laplacianos y Geometría Espectral](#22-operadores-laplacianos-y-geometría-espectral)  
   2.3 [Embedding Espectral y Reducción de Dimensionalidad](#23-embedding-espectral-y-reducción-de-dimensionalidad)  
   2.4 [Estabilidad Espectral](#24-estabilidad-espectral)  
3. [Metodología y Análisis Estático](#3-metodología-y-análisis-estático)

   3.1 [Selección de Datos y Configuración de Parámetros](#31-selección-de-datos-y-configuración-de-parámetros)  
   3.2 [Análisis de la Estructura Estática (Benchmark)](#32-análisis-de-la-estructura-estática-benchmark)  
   3.3 [Métricas de Validación](#33-métricas-de-validación)  
4. [Análisis Dinámico y Transiciones de Régimen](#4-análisis-dinámico-y-transiciones-de-régimen)

   4.1 [Protocolo de Ventanas Móviles (Rolling Windows)](#41-protocolo-de-ventanas-móviles-rolling-windows)  
   4.2 [Estabilidad Temporal y Alineación de Procrustes](#42-estabilidad-temporal-y-alineación-de-procrustes)  
   4.3 [Análisis de la Evolución Espectral](#43-análisis-de-la-evolución-espectral)  
   4.4 [Identificación de Regímenes](#44-identificación-de-regímenes)  
5. [Conclusiones](#5-conclusiones)

---

## 1. Introducción
Contexto del mercado como sistema complejo y el objetivo de identificar transiciones de régimen mediante geometría espectral.

## 2. Marco Teórico

### 2.1 Mercado financiero como red compleja
Ya redactado.

### 2.2 Operadores Laplacianos y Geometría Espectral
Definición de $L_{sym}$ y justificación de su uso para capturar la topología. Interpretación del espectro.

### 2.3 Embedding Espectral y Reducción de Dimensionalidad
Fundamento de los autovectores ($u_2, u_3$) para preservar proximidad local (Belkin & Niyogi).

### 2.4 Estabilidad Espectral
El Teorema de Davis-Kahan como garantía de que los cambios en el embedding son estructurales y no ruido estadístico.

## 3. Metodología y Análisis Estático

### 3.1 Selección de Datos y Configuración de Parámetros
Descripción del universo de activos y justificación del parámetro $k$ en el kernel self-tuning.

### 3.2 Análisis de la Estructura Estática (Benchmark)
Aplicación de la teoría a una ventana fija para validar que el embedding recupera la estructura sectorial.

### 3.3 Métricas de Validación
Definición de la Información Mutua Normalizada (NMI) y el Gap Espectral como descriptores de la modularidad observada.

## 4. Análisis Dinámico y Transiciones de Régimen

### 4.1 Protocolo de Ventanas Móviles (Rolling Windows)
Definición de *window* y *step*.

### 4.2 Estabilidad Temporal y Alineación de Procrustes
Resolución de la ambigüedad ortogonal para permitir la comparación de embeddings sucesivos.

### 4.3 Análisis de la Evolución Espectral
Resultados de la fluctuación del Gap y NMI a lo largo del tiempo.

### 4.4 Identificación de Regímenes
Discusión sobre los estados de modularidad sectorial vs. acoplamiento global (énfasis en el periodo COVID-19).

## 5. Conclusiones
Síntesis de cómo el espectro del Laplaciano actúa como un termómetro de la cohesión del mercado.
