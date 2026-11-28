# Source de `Doc/presentation_parties_A_B.pptx`

Deck des Parties A et B sur le template `doc_ref/template-CentraleSupélec-PPT.pptx`
(python-pptx ; animations et transitions écrites en XML dans `anim.py`). Depuis ce dossier :

    mkdir layers img
    python -c "import zipfile; zipfile.ZipFile('../../../doc_ref/template-CentraleSupélec-PPT.pptx').extract('ppt/media/image2.png', 'tpl')"
    python extract_layers.py tpl/ppt/media/image2.png layers    # anneaux du logo découpés en calques tournants
    python thom_schema.py img/thom_schema.png
    python corner_schema.py img/corner_schema.png
    python build_deck.py ../../../doc_ref/template-CentraleSupélec-PPT.pptx ../../../figures ../../presentation_parties_A_B.pptx

Paquets : python-pptx, lxml, pillow, numpy, scipy, matplotlib. Les figures viennent des notebooks A et B.
