
add_input_file_links()
dune_symlink_to_source_files(FILES injections.dat)


# compile MICP simplified chemistry column setup
dune_add_test(NAME biomineralization
              SOURCES main.cc
              COMPILE_DEFINITIONS TYPETAG=MICPColumnSimpleChemistryBoxTypeTag)
