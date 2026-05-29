# ESM-2 Pangénome Embeddings Explorer

Ce dépôt contient une interface de consultation pour explorer et comparer plusieurs méthodes de similarité entre protéines en utilisant des embeddings ESM-2. Le jeu de données et l'interface web sont pré-configurés et prêts à l'emploi. Il n'est pas conçu pour reconstruire la pipeline (lourde) de création des données.

## Installation et démarrage

### 1. Cloner le dépôt

```bash
git clone https://github.com/HugoBrtn/Etude_ESM-2.git
cd Etude_ESM-2
git lfs install
git lfs pull
```

### 2. Créer l'environnement Conda

**Windows (PowerShell)** :
```powershell
conda env create -f environment.yml
conda activate etude_esm2
```

**Linux / macOS** :
```bash
conda env create -f environment.yml
conda activate etude_esm2
```

### 3. Régénérer le bundle (optionnel)

Si le dépôt a été déplacé ou les données modifiées, régénérez le bundle :
```bash
python scripts/prepare_interface_data.py
```

 **Cette étape peut prendre 1-5 minutes** selon votre système. 

### 4. Lancer l'interface

```bash
python scripts/run_ui.py
```

Ouvrez votre navigateur sur :
```
http://127.0.0.1:5000
```

---


## Contexte scientifique

Ce projet fait partie d'une **étude sur la construction de multi-graphes de pangénomes**. L'objectif est d'identifier les meilleures post processing sur les embeddings ESM-2 pour capturer des similarités biologiques entre familles de gènes à travers différentes espèces.

Les approches classiques de similarité (séquence, structure) peuvent manquer des relations biologiques complexes ou distantes, particulièrement quand la conservation est faible. Les **embeddings ESM-2** permettent de capturer des similarités **fonctionnelles et structurales latentes**, souvent invisibles aux méthodes d'alignement traditionnelles (MMseqs2, Needleman-Wunsch, TM-align, Foldseek).

### Données utilisées

L'interface contient des données de protéines originaires de :
- ***Escherichia coli* (souche K-12 MG1655)**
- ***Bacillus subtilis* (souche 168)**

**Critères de sélection** :
- Structures disponible sur **AlphaFold DB** avec pLDDT moyen > 85.
- Couverture query et target > 0.3 dans les alignements MMseqs2.

### Méthodes de similarité comparées

L'interface compare quatre catégories de similarité :

#### 1. **Alignements de séquence**
- **MMseqs2** : similarité de séquence locale.
- **Needleman-Wunsch** : similarité de séquence globale.

#### 2. **Alignements structuraux**
- **TM-align** : score d'alignement structural basé sur superposition 3D.
- **Foldseek** : recherche ultrapide de similarité structurale.

#### 3. **Embeddings ESM-2 (30 méthodes de post-processing)**

Cinq **méthodes de pooling** testées pour agréger les embeddings de séquences :

| Pooling | Description |
|---------|-------------|
| `max` | Valeur maximale par dimension |
| `mean` | Moyenne des dimensions |
| `sum` | Somme des dimensions |
| `bos` | Embedding du token CLS (Beginning Of Sequence) |
| `mahalanobis` | Distance de Mahalanobis—sensible aux variations structurées |

Avec six **conditions post-traitement** pour gérer les dimensions outliers (valeurs extrêmes influençant les mesures) :

| Condition | Approche |
|-----------|----------|
| `raw` | Pas de traitement |
| `normalized` | Normalisation standard par dimension |
| `mean_outliers_filtered` | Dimensions outliers (sur moyennes) retirées |
| `mean_outliers_only` | Comparaison uniquement sur dimensions outliers des moyennes |
| `std_outliers_filtered` | Dimensions outliers (sur écarts-types) retirées |
| `std_outliers_only` | Comparaison uniquement sur dimensions outliers des écarts-types |

**Total : 5 pooling × 6 conditions = 30 variantes d'embeddings**

Détection des outliers : médiane + MAD (Median Absolute Deviation) ou IQR selon la distribution avec seuil à trois dimensiosn outliers.

**Mesure de similarité** : similarité cosinus sur les vecteurs d'embeddings normalisés.

---


## Utilisation de l'interface

L'interface permet de :
- **Rechercher des protéines** par accession ou nom
- **Afficher les paires** et leurs métriques de similarité (séquence, structure, embeddings)
- **Filtre par méthode** : comparhez MMseqs2 vs. embeddings vs. TM-align
- **Visualiser les alignements** : accédez aux fichiers d'alignement bruts
- **Consulter les structures** : affichage 3D des prédictions AlphaFold
- **Télécharger les données** : exportez les CSV globales

---

## Organisation du dépôt

```
.
├── README.md                      # Ce fichier
├── environment.yml                # Environnement Conda (Flask, numpy, scipy, statsmodels)
├── scripts/
│   ├── run_ui.py                  # Point d'entrée : lance le serveur Flask
│   ├── prepare_interface_data.py   # Régénère les bundles après modification de data/
│   └── global_aggregator.py        # Fonctions internes d'agrégation CSV
├── ui/
│   ├── app.py                      # Serveur Flask + endpoints API
│   ├── index.html                  # Interface HTML/CSS/JS
│   └── app_data.js                 # Bundle de données pré-calculé (Git LFS)
└── data/
    ├── GLOBAL/
    │   ├── proteins_global.csv     # Métadonnées par protéine
    │   └── pairs_global.csv        # Tableau des paires + métriques
    ├── inputs/                     # Protéines brutes (séquences, structures, embeddings)
    ├── alignment_mmseq2/           # Alignements MMseqs2
    ├── alignment_needleman_wunsh/  # Alignements Needleman-Wunsch
    ├── alignment_structure_tmscore/# Alignements TM-align
    ├── alignment_structure_foldseek/# Alignements Foldseek
    └── embedding_similarity/       # Similarités ESM-2 par pooling/condition
```
