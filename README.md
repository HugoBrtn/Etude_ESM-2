# comparaison_seq_struc_emb_new

Pipeline de comparaison de protéines pour séquences, structures et embeddings, avec une interface web dédiée.

## Objectif

La pipeline récupère des protéines depuis UniProt, vérifie si une structure est disponible dans AlphaFold DB, peut recalculer une structure avec ColabFold si nécessaire, calcule les embeddings ESM-2, produit les alignements de séquence et de structure, puis agrège tout dans des CSV globaux et dans l'UI.

## Vue D'ensemble

Flux principal:
1. Collecte UniProt + structure AlphaFold (ou ColabFold si besoin)
2. Embeddings ESM-2 multi-pooling
3. Alignements séquence (MMseqs2, puis Needleman-Wunsch filtré)
4. Alignements structure (TM-align, filtré)
5. Similarités embeddings (cosinus, filtré)
6. Agrégation CSV globaux
7. Bundle de données pour l'interface web

Comportement de reprise:
- Les scripts de calcul sautent une paire si le fichier final existe déjà.
- Le mode full proteome ajoute des checkpoints par étape (reprise globale).

Filtrage MMseqs2 (coverage):
- Deux seuils sont disponibles: un critère AND (qcov et tcov) et un critère OR (qcov ou tcov).
- Le mode peut appliquer aucun, un seul, ou les deux critères.
- Par défaut, les deux critères sont activés avec les valeurs de pipeline_config.env.

## Prérequis

- Environnement conda: `comparaison_emb`
- Outils: `mmseqs`, `TMalign`, `colabfold_batch`
- Paquets Python: `torch`, `esm`, `numpy`, `scipy`, `biopython`, `flask`

Créer ou réparer l'environnement complet:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/setup_comparaison_emb.sh
```

Vérifier l'environnement:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/check_env.sh
```

Réinstaller seulement les outils externes si besoin:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/install_tools.sh
```

## Mise En Route

Avant de lancer quoi que ce soit, active l'environnement conda:

```bash
conda activate comparaison_emb
python -c "import torch; print(f'PyTorch: {torch.__version__}'); print(f'CUDA available: {torch.cuda.is_available()}')"
```

Si `CUDA available: True`, l'étape embeddings utilisera le GPU.

**Si CUDA n'est pas disponible**, réinstalle PyTorch avec le support CUDA:

```bash
# Détecte automatiquement ta version CUDA et l'installe correctement
bash Code/comparaison_seq_struc_emb_new/scripts/setup_comparaison_emb.sh --gpu --recreate

# Ou manuellement avec conda:
conda activate comparaison_emb
conda install -y -c pytorch -c nvidia pytorch torchvision torchaudio pytorch-cuda=12.4
```

Puis vérifie à nouveau:

```bash
python -c "import torch; print(torch.cuda.is_available())"  # Doit afficher True
python -c "print(torch.cuda.get_device_name(0))"  # Affiche le GPU utilisé
```

## Étapes De La Pipeline

1. `collect` - collecte des protéines depuis UniProt / AlphaFold DB, filtre pLDDT, option ColabFold
2. `collect` - même logique pour la deuxième espèce
3. `embeddings` - embeddings ESM-2 multi-pooling
4. `mmseq2` - alignements MMseqs2
5. `nw` - alignements Needleman-Wunsch
6. `tm` - alignements structurels TM-align
7. `similarities` - similarités cosinus d'embeddings
8. `csv` - construction des CSV globaux
9. `ui` - génération du bundle de données pour l'interface

## Filtrage MMseqs2 (qcov/tcov)

Le filtrage s'applique aux étapes `nw`, `tm`, `similarities` et à l'agrégation globale.

Modes:
- `none`: pas de filtrage
- `and`: qcov ET tcov >= seuil AND
- `or`: qcov OU tcov >= seuil OR
- `both`: applique AND et OR (par défaut)

Paramètres (fichier: pipeline_config.env):
- `MMSEQ2_COVERAGE_AND_THRESHOLD`
- `MMSEQ2_COVERAGE_OR_THRESHOLD`
- `MMSEQ2_COVERAGE_MODE`

Exemple en ligne de commande:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_full_proteome.sh \
	--steps nw,tm,similarities \
	--mmseq2-coverage-and 0.2 \
	--mmseq2-coverage-or 0.6 \
	--mmseq2-coverage-mode both
```

## Guide De Réutilisation (Pipeline Personnalisée)

1. Préparer l'environnement (outils + conda).
2. Choisir un mode de collecte:
	 - taxon: une espèce via `--taxon-id`
	 - accessions explicites: liste d'IDs UniProt
3. Lancer les étapes dans l'ordre (collect -> embeddings -> alignements -> similarities -> csv -> ui).
4. Ouvrir l'UI si besoin.

Exemple minimal (taxon, 50 protéines):

```bash
conda activate comparaison_emb
bash Code/comparaison_seq_struc_emb_new/scripts/run_pipeline.sh \
	--taxon-id 83333 \
	--count 50
```

Recalcul ciblé d'étapes sans casser l'existant:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_full_proteome.sh \
	--steps embeddings,mmseq2 --force
```

## Script De Test

Le test lance un sous-ensemble de protéines par espèce (par défaut 20). En foreground:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_test_subset.sh
```

Définir le nombre de protéines:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_test_subset.sh --count 50
```

Définir le filtrage MMseqs2:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_test_subset.sh \
	--mmseq2-coverage-and 0.2 \
	--mmseq2-coverage-or 0.6 \
	--mmseq2-coverage-mode both
```

En arrière-plan (recommandé):

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_test_subset.sh --background
```

Le script affichera automatiquement les commandes de monitoring:

```bash
# Suivre la progression
tail -f Code/comparaison_seq_struc_emb_new/data/.logs/run_test_subset.nohup.log

# Voir l'utilisation GPU (met à jour toutes les 2 secondes)
watch -n 2 'nvidia-smi | grep -A 20 Processes'

# Arrêter le process
kill $(cat Code/comparaison_seq_struc_emb_new/data/.logs/run_test_subset.pid)
```

## Pipeline Full Proteome

Lancer toute la pipeline (foreground):

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_full_proteome.sh --full
```

Ou en arrière-plan:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_full_proteome.sh --full --background
```

Lancer seulement certaines étapes:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_full_proteome.sh \
	--steps collect,embeddings,mmseq2
```

Reprendre plus tard à partir de l'étape 5:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_full_proteome.sh \
	--steps nw,tm,similarities,csv,ui --background
```

Forcer le recalcul:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_full_proteome.sh \
	--steps embeddings,mmseq2 --force --background
```

Recalculer une structure manquante avec ColabFold:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_full_proteome.sh \
	--steps collect --colabfold
```

Ajuster le seuil pLDDT et les threads MMseq2:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_full_proteome.sh \
	--steps collect,embeddings,mmseq2 \
	--plddt-threshold 85 \
	--threads-mmseq2 16
```

Changer le filtrage MMseqs2:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_full_proteome.sh \
	--steps nw,tm,similarities \
	--mmseq2-coverage-and 0.2 \
	--mmseq2-coverage-or 0.6 \
	--mmseq2-coverage-mode both
```

## Parametres MMseqs2 Retenus

Les alignements sont lances avec ces options par defaut:

- `--prefilter-mode 1`: prefiltrage plus rapide pour reduire le temps de recherche.
- `--alignment-mode 3`: alignement local avec sortie detaillee.
- `--cov-mode 0`: couverture calculee sur la sequence complete.
- `--min-seq-id 0.1`: identite minimale (10%) pour ne pas jeter trop d'alignements faibles.
- `-e 10`: seuil d'e-value pour filtrer les alignements peu significatifs.
- `--max-seqs 25`: limite le nombre d'alignements conserves par requete.
- `--gap-open 11`: penalite d'ouverture de gap.
- `--gap-extend 1`: penalite d'extension de gap.
- `--format-mode 4`: format tabulaire compact pour parser facilement.

En plus, le filtrage MMseqs2 (qcov/tcov) est applique aux etapes `nw`, `tm`, `similarities` et a l'aggregation globale via:

- `MMSEQ2_COVERAGE_MODE`: mode de filtrage (`none`, `and`, `or`, `both`).
- `MMSEQ2_COVERAGE_AND_THRESHOLD`: seuil pour appliquer qcov ET tcov.
- `MMSEQ2_COVERAGE_OR_THRESHOLD`: seuil pour appliquer qcov OU tcov.

Surcharger depuis les launchers:

- `run_test_subset.sh` et `run_full_proteome.sh` acceptent `--mmseq2-coverage-and`, `--mmseq2-coverage-or`, `--mmseq2-coverage-mode`.
- `run_full_proteome.sh` accepte `--threads-mmseq2` pour ajuster les threads.
- Pour modifier d'autres options MMseqs2 (ex: `--min-seq-id`, `-e`, `--max-seqs`), il faut editer [Code/comparaison_seq_struc_emb_new/scripts/mmseqs2_align.py](Code/comparaison_seq_struc_emb_new/scripts/mmseqs2_align.py).

Noms d'étapes acceptés: `collect`, `embeddings`, `mmseq2`, `nw`, `tm`, `similarities`, `csv`, `ui`, `full`.

## Pipeline Personnalisée

20 protéines pour un taxon donné:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_pipeline.sh \
	--taxon-id 83333 \
	--count 20
```

Toutes les protéines d'un taxon:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_pipeline.sh \
	--taxon-id 83333 \
	--all \
	--colabfold
```

Lancer en arrière-plan:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_pipeline.sh \
	--taxon-id 83333 \
	--all \
	--colabfold \
	--background
```

## Comportement De La Collecte

La collecte tente d'abord de récupérer une structure depuis AlphaFold DB. Si la structure est absente et que `--colabfold-if-missing` est activé, ColabFold est utilisé pour recalculer une structure locale.

Le filtre structurel s'applique sur le pLDDT moyen extrait du PDB. Si `mean_plddt <= plddt_threshold`, la protéine est rejetée.

Exemple:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_full_proteome.sh \
	--steps collect \
	--colabfold \
	--plddt-threshold 85
```

## Exécution En Arrière-Plan

Tous les scripts launchers supportent `--background` pour une exécution détachée:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_test_subset.sh --background
bash Code/comparaison_seq_struc_emb_new/scripts/run_full_proteome.sh --full --background
bash Code/comparaison_seq_struc_emb_new/scripts/run_pipeline.sh --taxon-id 83333 --all --background
```

Avec `--background`, le script affiche automatiquement:
- **PID**: identificateur du process
- **Fichier PID**: chemin pour retrouver le PID plus tard
- **Commandes de monitoring**: prêtes à copier-coller

Exemple de sortie:

```
[INFO] Process started in background (PID: 12345)
[INFO] Main log: Code/comparaison_seq_struc_emb_new/data/.logs/run_test_subset.nohup.log
[INFO] Monitor progress:
       tail -f Code/comparaison_seq_struc_emb_new/data/.logs/run_test_subset.nohup.log
[INFO] Check GPU usage:
       watch -n 2 'nvidia-smi | grep -A 20 Processes'
[INFO] Stop pipeline:
       kill 12345  # or: kill $(cat Code/comparaison_seq_struc_emb_new/data/.logs/run_test_subset.pid)
[INFO] Check final status:
       ps -p 12345 && echo Running || echo Done
```

## Reprise Et Checkpoints

Le script full proteome crée des marqueurs dans `Code/comparaison_seq_struc_emb_new/data/.logs/.checkpoints`.

- Si une étape est déjà terminée, elle est sautée lors d'un relancement.
- `--force` permet de recalculer une étape même si un marqueur existe.
- Les logs d'étapes contiennent le temps de début, la durée de l'étape et le cumul.

## Données Produites

```text
data/
	inputs/
		<species_key>/<accession>/
			sequence.fasta
			structure.pdb
			structure_colabfold.pdb
			esm2_multipooling.pt
			metadata.json
	alignment_mmseq2/
	alignment_needleman_wunsh/
	alignment_structure_foldseek/
	embedding_similarity/
	GLOBAL/
		proteins_global.csv
		pairs_global.csv
```

## Fichiers Et Rôles (Scripts)

| Fichier | Rôle | Entrées / Sorties |
| --- | --- | --- |
| [Code/comparaison_seq_struc_emb_new/scripts/run_test_subset.sh](Code/comparaison_seq_struc_emb_new/scripts/run_test_subset.sh) | Pipeline de test (subset) | Lance toutes les étapes sur un petit échantillon |
| [Code/comparaison_seq_struc_emb_new/scripts/run_full_proteome.sh](Code/comparaison_seq_struc_emb_new/scripts/run_full_proteome.sh) | Pipeline full proteome | Orchestration par étapes + checkpoints |
| [Code/comparaison_seq_struc_emb_new/scripts/run_pipeline.sh](Code/comparaison_seq_struc_emb_new/scripts/run_pipeline.sh) | Pipeline par taxon | Collecte taxon + étapes aval |
| [Code/comparaison_seq_struc_emb_new/scripts/collect.py](Code/comparaison_seq_struc_emb_new/scripts/collect.py) | Collecte UniProt/AlphaFold/ColabFold | Sorties: FASTA, PDB, metadata |
| [Code/comparaison_seq_struc_emb_new/scripts/compute_esm2_embeddings.py](Code/comparaison_seq_struc_emb_new/scripts/compute_esm2_embeddings.py) | Embeddings ESM-2 multi-pooling | Sortie: esm2_multipooling.pt |
| [Code/comparaison_seq_struc_emb_new/scripts/mmseqs2_align.py](Code/comparaison_seq_struc_emb_new/scripts/mmseqs2_align.py) | Alignements MMseqs2 (DB-first, index cache, fallback pairwise) | Sorties: pairwise_summary.csv, alignment.txt.gz |
| [Code/comparaison_seq_struc_emb_new/scripts/needleman_align.py](Code/comparaison_seq_struc_emb_new/scripts/needleman_align.py) | Alignements Needleman-Wunsch | Sorties: pairwise_summary.csv, alignment.txt.gz |
| [Code/comparaison_seq_struc_emb_new/scripts/compare_structure_tmalign.py](Code/comparaison_seq_struc_emb_new/scripts/compare_structure_tmalign.py) | Alignements structure TM-align | Sorties: pairwise_summary.csv, alignment.txt.gz |
| [Code/comparaison_seq_struc_emb_new/scripts/compute_embedding_similarity.py](Code/comparaison_seq_struc_emb_new/scripts/compute_embedding_similarity.py) | Similarites embeddings | Sortie: pairwise_summary.csv |
| [Code/comparaison_seq_struc_emb_new/scripts/global_aggregator.py](Code/comparaison_seq_struc_emb_new/scripts/global_aggregator.py) | Agrégation CSV + bundle UI | Sorties: proteins_global.csv, pairs_global.csv, app_data.js |
| [Code/comparaison_seq_struc_emb_new/scripts/build_global_files.py](Code/comparaison_seq_struc_emb_new/scripts/build_global_files.py) | Rebuild CSV globaux | Sorties: proteins_global.csv, pairs_global.csv |
| [Code/comparaison_seq_struc_emb_new/scripts/build_ui_data.py](Code/comparaison_seq_struc_emb_new/scripts/build_ui_data.py) | Rebuild bundle UI | Sortie: app_data.js |
| [Code/comparaison_seq_struc_emb_new/scripts/pipeline_config.env](Code/comparaison_seq_struc_emb_new/scripts/pipeline_config.env) | Seuils globaux | Variables d'environnement (ex: couverture MMseqs2) |
| [Code/comparaison_seq_struc_emb_new/scripts/check_env.sh](Code/comparaison_seq_struc_emb_new/scripts/check_env.sh) | Diagnostic env | Vérifie outils et Python |
| [Code/comparaison_seq_struc_emb_new/scripts/setup_comparaison_emb.sh](Code/comparaison_seq_struc_emb_new/scripts/setup_comparaison_emb.sh) | Setup conda | Installe conda + deps |
| [Code/comparaison_seq_struc_emb_new/scripts/install_tools.sh](Code/comparaison_seq_struc_emb_new/scripts/install_tools.sh) | Install outils externes | MMseqs2, TMalign, Foldseek, etc. |
| [Code/comparaison_seq_struc_emb_new/scripts/diagnose_cuda.py](Code/comparaison_seq_struc_emb_new/scripts/diagnose_cuda.py) | Diagnostic GPU | Affiche CUDA/Torch |
| [Code/comparaison_seq_struc_emb_new/scripts/foldseek_structure.py](Code/comparaison_seq_struc_emb_new/scripts/foldseek_structure.py) | Alignements structure Foldseek (forward + reverse, robust handling) | Optionnel, utilisé si Foldseek est installé |
| [Code/comparaison_seq_struc_emb_new/scripts/cleanup_storage.py](Code/comparaison_seq_struc_emb_new/scripts/cleanup_storage.py) | Nettoyage des fichiers temporaires et stockage | Supprime anciens logs / cache |

## Fichiers Et Rôles (UI)

| Fichier | Rôle |
| --- | --- |
| [Code/comparaison_seq_struc_emb_new/ui/app.py](Code/comparaison_seq_struc_emb_new/ui/app.py) | Serveur Flask de l'interface |
| [Code/comparaison_seq_struc_emb_new/ui/index.html](Code/comparaison_seq_struc_emb_new/ui/index.html) | Frontend statique |
| [Code/comparaison_seq_struc_emb_new/ui/app_data.js](Code/comparaison_seq_struc_emb_new/ui/app_data.js) | Données injectées (générées) |
```

## Fichiers Globaux

- `proteins_global.csv`: métadonnées par protéine, annotations, type de sous-unité, chemins fasta/pdb
- `pairs_global.csv`: métriques par paire, incluant MMseqs2, Needleman, Foldseek/TM-like et embeddings

## Interface Web

Lancer l'UI:

```bash
python Code/comparaison_seq_struc_emb_new/ui/app.py
```

Ouvrir:

```text
http://localhost:5000
```

L'interface affiche les paires, les annotations, la FASTA, le visualiseur 3D et les alignements.

## Conseils De Performance

- Utiliser un GPU pour `compute_esm2_embeddings.py` accélère fortement la pipeline.
- L'étape embeddings utilise maintenant le batching, `torch.inference_mode()` et `autocast` sur CUDA.
- `MMseqs2` et `Foldseek` utilisent par défaut tous les cœurs disponibles.
- Les scripts de lancement évitent les lectures sur stdin quand ils sont lancés via `nohup`.

## Résumé Rapide

```bash
conda activate comparaison_emb
bash Code/comparaison_seq_struc_emb_new/scripts/run_full_proteome.sh --full
```

Ou pour un test rapide:

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/run_test_subset.sh
```

## Changements récents (Mai 2026)

Les évolutions récentes apportent des options de filtrage par relation d'espèces et l'affichage du coefficient de Spearman dans l'interface :

- UI: nouveau contrôle `Relation espèces` (valeurs: `Les deux`, `Intra-espèces`, `Inter-espèces`). Permet d'afficher uniquement les paires intra-espèces (même espèce), inter-espèces (espèces différentes) ou les deux.
- Backend UI/API: l'endpoint LOESS renvoie maintenant `pearson` et `spearman` (utilisé pour afficher `r` et `ρ` dans l'en-tête du plot).
- Scripts: ajout d'un argument CLI `--pair-species-mode` desservant les scripts de comparaison par paire suivants :
	- `scripts/mmseqs2_align.py`
	- `scripts/needleman_align.py`
	- `scripts/foldseek_structure.py`
	- `scripts/compare_structure_tmalign.py`
	- `scripts/compute_embedding_similarity.py`

	Valeurs acceptées: `both` (défaut), `intra`, `inter`.

- Global aggregator: `pairs_global.csv` contient déjà `query_species` et `target_species` qui servent le filtrage côté UI.

Utilisation CLI exemple (calculer uniquement paires intra-espèces pour MMseqs2):

```bash
bash Code/comparaison_seq_struc_emb_new/scripts/mmseqs2_align.py --pair-species-mode intra
```

Notes:
- Par défaut (`both`), le comportement existant est conservé.
- L'option peut être passée via la variable d'environnement `PAIR_SPECIES_MODE` si vous préférez exporter l'option globalement avant d'appeler plusieurs scripts.

Remarques récentes:
- Nouveaux flags CLI disponibles: `--ultra-fast`, `--use-gpu`, et bascule DB pour MMseqs2 via `--mmseqs-use-db` / `--mmseqs-no-db` (exposés depuis les wrappers `run_test_subset.sh` et `run_full_proteome.sh`).
- Foldseek: intégré dans les scripts et ajouté aux helpers d'installation; le wrapper gère désormais les runs "reverse" manquants sans échouer la pipeline.
- UI / agrégation: le champ `bits` produit par certains wrappers est en cours de retrait de l'UI — l'agrégateur peut encore contenir des colonnes `mmseq2_bits` / `foldseek_bits` (valeurs vides/NA) ; une passe de nettoyage est recommandée avant un déploiement final.

