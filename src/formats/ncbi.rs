use std::collections::HashMap;
use std::fmt::{Debug, Display};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::str::FromStr;

use crate::base::{GeneralTaxonomy, NodeMetadata, TaxonomyName};
use crate::errors::{Error, ErrorKind, TaxonomyResult};
use crate::rank::TaxRank;
use crate::taxonomy::Taxonomy;

const NODES_FILENAME: &str = "nodes.dmp";
const NAMES_FILENAME: &str = "names.dmp";

/// Loads a NCBI taxonomy from the given directory.
/// The directory should contain at least two files: `nodes.dmp` and `names.dmp`.
pub fn load<P: AsRef<Path>>(ncbi_directory: P) -> TaxonomyResult<GeneralTaxonomy> {
    let dir = ncbi_directory.as_ref();
    let nodes_file = std::fs::File::open(dir.join(NODES_FILENAME))?;
    let names_file = std::fs::File::open(dir.join(NAMES_FILENAME))?;

    // First we go through the nodes
    let mut tax_ids: Vec<String> = Vec::new();
    let mut parents: Vec<String> = Vec::new();
    let mut ranks: Vec<TaxRank> = Vec::new();
    let mut node_metadata: Vec<NodeMetadata> = Vec::new();
    let mut tax_to_idx: HashMap<String, usize> = HashMap::new();

    for (ix, line) in BufReader::new(nodes_file).lines().enumerate() {
        let line_str = line?;
        let mut fields: Vec<String> = line_str.split("\t|\t").map(|x| x.to_string()).collect();

        // Clean up the last field which has a trailing "\t|" or "\t|\n"
        if let Some(last) = fields.last_mut() {
            *last = last.trim_end_matches("\t|").trim().to_string();
        }

        if fields.len() < 13 {
            // Should be at least 13 fields (original format), possibly up to 18 for newer dumps
            return Err(Error::new(ErrorKind::ImportError {
                line: ix,
                msg: format!("Not enough fields in nodes.dmp (got {}); bad line?", fields.len()),
            }));
        }

        let tax_id = fields[0].trim().to_string();
        let parent_tax_id = fields[1].trim().to_string();
        let rank = fields[2].trim().to_string();

        // Parse node metadata (fields 3-17)
        let metadata = NodeMetadata {
            embl_code: fields.get(3).map(|s| s.trim().to_string()).unwrap_or_default(),
            division_id: fields.get(4).map(|s| s.trim().to_string()).unwrap_or_default(),
            inherited_div_flag: fields.get(5).map(|s| s.trim() == "1").unwrap_or(false),
            genetic_code_id: fields.get(6).map(|s| s.trim().to_string()).unwrap_or_default(),
            inherited_gc_flag: fields.get(7).map(|s| s.trim() == "1").unwrap_or(false),
            mitochondrial_genetic_code_id: fields.get(8).map(|s| s.trim().to_string()).unwrap_or_default(),
            inherited_mgc_flag: fields.get(9).map(|s| s.trim() == "1").unwrap_or(false),
            genbank_hidden_flag: fields.get(10).map(|s| s.trim() == "1").unwrap_or(false),
            hidden_subtree_root_flag: fields.get(11).map(|s| s.trim() == "1").unwrap_or(false),
            comments: fields.get(12).map(|s| s.trim().to_string()).unwrap_or_default(),
            plastid_genetic_code_id: fields.get(13).map(|s| s.trim().to_string()).unwrap_or_default(),
            inherited_pgc_flag: fields.get(14).map(|s| s.trim() == "1").unwrap_or(false),
            specified_species: fields.get(15).map(|s| s.trim().to_string()).unwrap_or_default(),
            hydrogenosome_genetic_code_id: fields.get(16).map(|s| s.trim().to_string()).unwrap_or_default(),
            inherited_hgc_flag: fields.get(17).map(|s| s.trim() == "1").unwrap_or(false),
        };

        tax_ids.push(tax_id.clone());
        parents.push(parent_tax_id);
        ranks.push(TaxRank::from_str(&rank)?);
        node_metadata.push(metadata);
        tax_to_idx.insert(tax_id, ix);
    }

    // TODO: fixme? this fails if we have unmapped parent nodes (i.e. file is truncated?)
    let mut parent_ids = Vec::with_capacity(parents.len());
    for (i, parent) in parents.into_iter().enumerate() {
        if let Some(idx) = tax_to_idx.get(&parent) {
            parent_ids.push(*idx);
        } else {
            return Err(Error::new(ErrorKind::ImportError {
                line: i + 1,
                msg: format!("Parent ID {} could not be found in nodes.dmp", parent),
            }));
        }
    }

    // And then grab their names by their idx
    let mut names: Vec<String> = vec![String::new(); tax_ids.len()];
    let mut all_names: Vec<TaxonomyName> = Vec::new();

    for (ix, line) in BufReader::new(names_file).lines().enumerate() {
        let line_str = line?;
        let mut fields: Vec<String> = line_str.split("\t|\t").map(|x| x.to_string()).collect();

        // Clean up the last field which has a trailing "\t|" or "\t|\n"
        if let Some(last) = fields.last_mut() {
            *last = last.trim_end_matches("\t|").trim().to_string();
        }

        if fields.len() < 4 {
            return Err(Error::new(ErrorKind::ImportError {
                line: ix,
                msg: format!("Not enough fields in names.dmp (got {}); bad line?", fields.len()),
            }));
        }

        let tax_id = fields[0].trim().to_string();
        let name = fields[1].trim().to_string();
        let unique_name = fields[2].trim().to_string();
        let name_class = fields[3].trim().to_string();

        let tax_id_idx = match tax_to_idx.get(&tax_id) {
            Some(&idx) => idx,
            None => {
                return Err(Error::new(ErrorKind::ImportError {
                    line: ix,
                    msg: format!("Tax ID {} in names.dmp not found in nodes.dmp", tax_id),
                }));
            }
        };

        if name_class.starts_with("scientific name") {
            // Store scientific name in the names vector
            names[tax_id_idx] = name;
        } else {
            // Store all other names in all_names
            all_names.push(TaxonomyName {
                tax_id_index: tax_id_idx,
                name,
                unique_name,
                name_class,
            });
        }
    }

    let mut gt =
        GeneralTaxonomy::from_arrays(tax_ids, parent_ids, Some(names), Some(ranks), None, None)?;

    // Set the additional NCBI-specific fields
    gt.all_names = all_names;
    gt.node_metadata = node_metadata;

    // Rebuild indices to include all_names in name_lookup
    gt.index();

    gt.validate_uniqueness()?;
    Ok(gt)
}

pub fn save<'t, T: 't, P: AsRef<Path>, X: Taxonomy<'t, T>>(
    tax: &'t X,
    out_dir: P,
) -> TaxonomyResult<()>
where
    T: Clone + Debug + Display + PartialEq,
{
    // This is the generic save function that works with any Taxonomy trait implementation
    // It doesn't have access to the extended NCBI fields, so it writes with default values
    let dir = out_dir.as_ref();
    std::fs::create_dir_all(dir)?;
    let mut node_writer = BufWriter::new(std::fs::File::create(dir.join(NODES_FILENAME))?);
    let mut name_writer = BufWriter::new(std::fs::File::create(dir.join(NAMES_FILENAME))?);

    let root = tax.root();
    for key in tax.traverse(root.clone())?.filter(|x| x.1).map(|x| x.0) {
        let name = tax.name(key.clone())?;
        let rank = tax.rank(key.clone())?;
        let parent = if key == root {
            format!("{}", key)
        } else {
            tax.parent(key.clone())?
                .map(|(x, _)| format!("{}", x))
                .unwrap_or_default()
        };

        // Write scientific name
        name_writer
            .write_all(format!("{}\t|\t{}\t|\t\t|\tscientific name\t|\n", &key, name).as_bytes())?;

        // Write node with all 18 fields (using defaults for extended fields)
        node_writer.write_all(
            format!(
                "{}\t|\t{}\t|\t{}\t|\t\t|\t\t|\t0\t|\t\t|\t0\t|\t\t|\t0\t|\t0\t|\t0\t|\t\t|\t\t|\t0\t|\t\t|\t\t|\t0\t|\n",
                &key,
                parent,
                rank.to_ncbi_rank(),
            )
            .as_bytes(),
        )?;
    }

    Ok(())
}

/// Save a GeneralTaxonomy to NCBI format with full metadata
/// This specialized version preserves all NCBI-specific fields
pub fn save_with_metadata(
    tax: &GeneralTaxonomy,
    out_dir: impl AsRef<Path>,
) -> TaxonomyResult<()> {
    let dir = out_dir.as_ref();
    std::fs::create_dir_all(dir)?;
    let mut node_writer = BufWriter::new(std::fs::File::create(dir.join(NODES_FILENAME))?);
    let mut name_writer = BufWriter::new(std::fs::File::create(dir.join(NAMES_FILENAME))?);

    // Write nodes.dmp with all metadata
    for (idx, tax_id) in tax.tax_ids.iter().enumerate() {
        let parent_id = if idx == 0 {
            tax_id.clone()
        } else {
            tax.tax_ids[tax.parent_ids[idx]].clone()
        };

        let rank = &tax.ranks[idx];
        let metadata = &tax.node_metadata[idx];

        // Convert boolean flags to "1" or "0"
        let to_flag = |b: bool| if b { "1" } else { "0" };

        node_writer.write_all(
            format!(
                "{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\t{}\t|\n",
                tax_id,
                parent_id,
                rank.to_ncbi_rank(),
                metadata.embl_code,
                metadata.division_id,
                to_flag(metadata.inherited_div_flag),
                metadata.genetic_code_id,
                to_flag(metadata.inherited_gc_flag),
                metadata.mitochondrial_genetic_code_id,
                to_flag(metadata.inherited_mgc_flag),
                to_flag(metadata.genbank_hidden_flag),
                to_flag(metadata.hidden_subtree_root_flag),
                metadata.comments,
                metadata.plastid_genetic_code_id,
                to_flag(metadata.inherited_pgc_flag),
                metadata.specified_species,
                metadata.hydrogenosome_genetic_code_id,
                to_flag(metadata.inherited_hgc_flag),
            )
            .as_bytes(),
        )?;

        // Write scientific name for this node
        let sci_name = &tax.names[idx];
        name_writer.write_all(
            format!("{}\t|\t{}\t|\t\t|\tscientific name\t|\n", tax_id, sci_name).as_bytes(),
        )?;
    }

    // Write all non-scientific names
    for tax_name in &tax.all_names {
        let tax_id = &tax.tax_ids[tax_name.tax_id_index];
        name_writer.write_all(
            format!(
                "{}\t|\t{}\t|\t{}\t|\t{}\t|\n",
                tax_id, tax_name.name, tax_name.unique_name, tax_name.name_class
            )
            .as_bytes(),
        )?;
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::taxonomy::Taxonomy;
    use std::io::Write;
    use tempfile::tempdir;

    #[test]
    fn can_import_ncbi() {
        let nodes =
            "1\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\t\t|\t0\t|\t\t|\t\t|\t0\t|
    10239\t|\t1\t|\tno_rank\t|\t\t|\t9\t|\t0\t|\t1\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\t\t|\t0\t|\t\t|\t\t|\t0\t|
    2\t|\t131567\t|\tsuperkingdom\t|\t\t|\t0\t|\t0\t|\t11\t|\t0\t|\t0\t|\t0\t|\t0\t|\t0\t|\t\t|\t\t|\t0\t|\t\t|\t\t|\t0\t|
    543\t|\t91347\t|\tfamily\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|\t\t|\t0\t|\t\t|\t\t|\t0\t|
    561\t|\t543\t|\tgenus\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|\t\t|\t0\t|\t\t|\t\t|\t0\t|
    562\t|\t561\t|\tspecies\t|\tEC\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t1\t|\t0\t|\t\t|\t\t|\t0\t|\t1\t|\t\t|\t0\t|
    1224\t|\t2\t|\tphylum\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|\t\t|\t0\t|\t\t|\t\t|\t0\t|
    1236\t|\t1224\t|\tclass\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|\t\t|\t0\t|\t\t|\t\t|\t0\t|
    91347\t|\t1236\t|\torder\t|\t\t|\t0\t|\t1\t|\t11\t|\t1\t|\t0\t|\t1\t|\t0\t|\t0\t|\t\t|\t\t|\t0\t|\t\t|\t\t|\t0\t|
    131567\t|\t1\t|\tno rank\t|\t\t|\t8\t|\t1\t|\t1\t|\t1\t|\t0\t|\t1\t|\t1\t|\t0\t|\t\t|\t\t|\t0\t|\t\t|\t\t|\t0\t|";
        let names = "1\t|\tall\t|\t\t|\tsynonym\t|
    1\t|\troot\t|\t\t|\tscientific name\t|
    10239\t|\tViruses\t|\t\t|\tscientific name\t|
    2\t|\tBacteria\t|\tBacteria <prokaryotes>\t|\tscientific name\t|
    543\t|\tEnterobacteriaceae\t|\t\t|\tscientific name\t|
    561\t|\tEscherchia\t|\t\t|\tmisspelling\t|
    561\t|\tEscherichia\t|\t\t|\tscientific name\t|
    561\t|\tEscherichia Castellani and Chalmers 1919\t|\t\t|\tauthority\t|
    562\t|\t\"Bacillus coli\" Migula 1895\t|\t\t|\tauthority\t|
    562\t|\t\"Bacterium coli commune\" Escherich 1885\t|\t\t|\tauthority\t|
    562\t|\t\"Bacterium coli\" (Migula 1895) Lehmann and Neumann 1896\t|\t\t|\tauthority\t|
    562\t|\tATCC 11775\t|\t\t|\ttype material\t|
    562\t|\tBacillus coli\t|\t\t|\tsynonym\t|
    562\t|\tBacterium coli\t|\t\t|\tsynonym\t|
    562\t|\tBacterium coli commune\t|\t\t|\tsynonym\t|
    562\t|\tCCUG 24\t|\t\t|\ttype material\t|
    562\t|\tCCUG 29300\t|\t\t|\ttype material\t|
    562\t|\tCIP 54.8\t|\t\t|\ttype material\t|
    562\t|\tDSM 30083\t|\t\t|\ttype material\t|
    562\t|\tE. coli\t|\t\t|\tcommon name\t|
    562\t|\tEnterococcus coli\t|\t\t|\tsynonym\t|
    562\t|\tEscherchia coli\t|\t\t|\tmisspelling\t|
    562\t|\tEscherichia coli\t|\t\t|\tscientific name\t|
    562\t|\tEscherichia coli (Migula 1895) Castellani and Chalmers 1919\t|\t\t|\tauthority\t|
    562\t|\tEscherichia sp. 3_2_53FAA\t|\t\t|\tincludes\t|
    562\t|\tEscherichia sp. MAR\t|\t\t|\tincludes\t|
    562\t|\tEscherichia/Shigella coli\t|\t\t|\tequivalent name\t|
    562\t|\tEschericia coli\t|\t\t|\tmisspelling\t|
    562\t|\tJCM 1649\t|\t\t|\ttype material\t|
    562\t|\tLMG 2092\t|\t\t|\ttype material\t|
    562\t|\tNBRC 102203\t|\t\t|\ttype material\t|
    562\t|\tNCCB 54008\t|\t\t|\ttype material\t|
    562\t|\tNCTC 9001\t|\t\t|\ttype material\t|
    562\t|\tbacterium 10a\t|\t\t|\tincludes\t|
    562\t|\tbacterium E3\t|\t\t|\tincludes\t|
    1224\t|\tAlphaproteobacteraeota\t|\t\t|\tsynonym\t|
    1224\t|\tAlphaproteobacteraeota Oren et al. 2015\t|\t\t|\tauthority\t|
    1224\t|\tProteobacteria\t|\t\t|\tscientific name\t|
    1224\t|\tProteobacteria Garrity et al. 2005\t|\t\t|\tauthority\t|
    1224\t|\tProteobacteria [class] Stackebrandt et al. 1988\t|\t\t|\tauthority\t|
    1224\t|\tnot Proteobacteria Cavalier-Smith 2002\t|\t\t|\tauthority\t|
    1224\t|\tproteobacteria\t|\tproteobacteria<blast1224>\t|\tblast name\t|
    1224\t|\tpurple bacteria\t|\t\t|\tcommon name\t|
    1224\t|\tpurple bacteria and relatives\t|\t\t|\tcommon name\t|
    1224\t|\tpurple non-sulfur bacteria\t|\t\t|\tcommon name\t|
    1224\t|\tpurple photosynthetic bacteria\t|\t\t|\tcommon name\t|
    1224\t|\tpurple photosynthetic bacteria and relatives\t|\t\t|\tcommon name\t|
    1236\t|\tGammaproteobacteria\t|\t\t|\tscientific name\t|
    91347\t|\tEnterobacterales\t|\t\t|\tscientific name\t|
    131567\t|\tbiota\t|\t\t|\tsynonym\t|
    131567\t|\tcellular organisms\t|\t\t|\tscientific name\t|";

        let dir = tempdir().unwrap();
        let path = dir.path();
        let mut nodes_file = std::fs::File::create(path.join(NODES_FILENAME)).unwrap();
        writeln!(nodes_file, "{}", nodes).unwrap();
        let mut names_file = std::fs::File::create(path.join(NAMES_FILENAME)).unwrap();
        writeln!(names_file, "{}", names).unwrap();

        let tax = load(path).unwrap();

        // Check basic taxonomy operations
        assert_eq!(
            Taxonomy::<&str>::name(&tax, "562").unwrap(),
            "Escherichia coli"
        );
        assert_eq!(
            Taxonomy::<&str>::rank(&tax, "562").unwrap(),
            TaxRank::Species
        );
        assert_eq!(
            Taxonomy::<&str>::children(&tax, "561").unwrap(),
            vec!["562"]
        );
        assert_eq!(
            Taxonomy::<&str>::parent(&tax, "562").unwrap(),
            Some(("561", 1.))
        );

        // Check that non-scientific names are loaded
        assert!(tax.all_names.len() > 0, "Should have loaded non-scientific names");

        // Check that we can find E. coli by its common name "E. coli"
        let found_by_common = tax.find_all_by_name("E. coli");
        assert!(found_by_common.contains(&"562"), "Should find E. coli by common name");

        // Check that we can find by synonym "Bacillus coli"
        let found_by_synonym = tax.find_all_by_name("Bacillus coli");
        assert!(found_by_synonym.contains(&"562"), "Should find E. coli by synonym");

        // Check node metadata is loaded (E. coli has embl_code "EC" and specified_species "1")
        let ecoli_idx = tax.to_internal_index("562").unwrap();
        assert_eq!(tax.node_metadata[ecoli_idx].embl_code, "EC");
        assert_eq!(tax.node_metadata[ecoli_idx].specified_species, "1");

        // Test save_with_metadata to preserve all fields
        let out = path.join("out");
        save_with_metadata(&tax, &out).unwrap();

        // now load again and validate a few taxids
        let tax2 = load(&out).unwrap();

        // Check E. coli (562)
        assert_eq!(
            Taxonomy::<&str>::name(&tax2, "562").unwrap(),
            "Escherichia coli"
        );
        assert_eq!(
            Taxonomy::<&str>::rank(&tax2, "562").unwrap(),
            TaxRank::Species
        );
        assert_eq!(
            Taxonomy::<&str>::parent(&tax2, "562").unwrap(),
            Some(("561", 1.))
        );

        // Check Escherichia (561)
        assert_eq!(Taxonomy::<&str>::name(&tax2, "561").unwrap(), "Escherichia");
        assert_eq!(
            Taxonomy::<&str>::rank(&tax2, "561").unwrap(),
            TaxRank::Genus
        );
        assert_eq!(
            Taxonomy::<&str>::parent(&tax2, "561").unwrap(),
            Some(("543", 1.))
        );

        // Check root (1)
        assert_eq!(Taxonomy::<&str>::name(&tax2, "1").unwrap(), "root");

        // Check children relationship preserved
        assert_eq!(
            Taxonomy::<&str>::children(&tax2, "561").unwrap(),
            vec!["562"]
        );

        // Check that non-scientific names survived round-trip
        assert!(tax2.all_names.len() > 0, "Non-scientific names should be preserved");
        let found_by_common2 = tax2.find_all_by_name("E. coli");
        assert!(found_by_common2.contains(&"562"), "Should still find E. coli by common name after round-trip");

        // Check metadata survived round-trip
        let ecoli_idx2 = tax2.to_internal_index("562").unwrap();
        assert_eq!(tax2.node_metadata[ecoli_idx2].embl_code, "EC", "EMBL code should be preserved");
        assert_eq!(tax2.node_metadata[ecoli_idx2].specified_species, "1", "Specified species flag should be preserved");
    }
}
