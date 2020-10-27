---
layout: post
date: 2018-04-24
title: "L'importation massive de riz au Bénin, une situation inquiétante"
categories: [Benin, Development]
tags: ['benin', 'opendata', 'riz', 'fr']
comments: true
highlight: true
toc: false
style: ["assets/css/posts/rice/trade.css", "assets/css/posts/general.css"]
---

### Un déficit commercial qui ne cesse de se creuser

Le Bénin a une balance commerciale déficitaire, et pas des moindres... Alors que la valeur monétaire des importations augmente exponentiellement chaque année, peu a changé en ce qui concerne les exportations. En 2016, le Bénin importait  \$2.63 milliards de biens, 6 fois plus que le total de ses exportations (\$410M).
  
<div class="vizdiv" id='tradeplot'></div>

À en croire ces données, la production locale est presque inexistante et la filière coton, dont les promesses de reviviscence nous ont été servies par les gourvenements successifs de ces dernières années, est `bien partie` pour boucher ce trou dans l'économie. 


Jusque là, je ne vous apprends rien de nouveau! Et puis, avoir un commerce déficitaire n'est pas forcément signe d'une économie mal en point<sup> [1](https://www.nytimes.com/2016/12/02/upshot/want-to-rev-up-the-economy-dont-worry-about-the-trade-deficit.html?smid=pl-share&_r=1), [2](https://www.contrepoints.org/2011/03/31/19383-deficit-de-la-balance-commerciale-et-alors)</sup>. Après tout, nous pourrions nous féliciter d'appartenir au [*cercle restreint* de pays](https://tradingeconomics.com/country-list/balance-of-trade) comme les États-Unis, la France, les UK, etc. Pourquoi s'alarmer donc ? *Right ?*


 
### Le riz, un problème qui n'inquiète personne.

**`$774M`**, c'est le montant total des importations de riz au Bénin en 2016, soit plus que les USA (\$714M), un pays avec une population 30 fois plus importante. Ce chiffre correspond à 1.5 millions de tonnes de riz et exclu les ré-exportations "officielles" vers d'autres pays. Comparé à l'importation de riz dans des pays avec une population similaire (Belgique, Sénégal, Zimbabwe, etc), cette quantité est inquiètante. 

<div class="vizdiv" id="trademap"> </div>

Il est vrai que nos habitudes alimentaires ont changé. Le riz, autrefois mets occasionnel, est devenu l'aliment de base de la population urbaine en expansion. Les raisons sont bien diverses, mais j'accuse entre autres, les programmes d'aide alimentaire distribuant leurs sacs de riz comme des petits pains. [Encore récemment](http://ortb.bj/index.php/societe/529-le-japon-fait-don-de-7000-tonnes-de-riz-au-benin), le Bénin a reçu du Japon, un *don de riz* pour renforcer la sécurité alimentaire. La commercialisation à bas prix de ce riz, au profit des couches défavorisées, pourraient générer un revenu utilisable pour réaliser des projets socio-communautaires. À première vue, il faudrait donc saluer l'initiative, mais considérant la **mauvaise gestion légendaire** de nos autorités étatiques, cette `coopération nippo-beninoise` n'aura pour seul résultat qu'alimenter notre dépendance au riz.

<div class="vizdiv" id="geodata"> </div>

Malgré l'importante consommation de riz blanc au Bénin, il serait quand même pernicieux d'insunier que le béninois consomme en moyenne 136 kg de riz par an (soit [trois plats de riz, par personne, chaque jour](https://health.gov/dietaryguidelines/2015/guidelines/appendix-3/)). Une analyse plus approfondie de la quantité de riz importée *per capita* dans les pays limitrophes, suggère qu'une partie serait écoulée par voies informelles (et donc non déclarée) au Togo et/ou au Nigeria. Étant donné qu'au Togo, [l'importation couvre les besoins non comblés par la production locale](https://riceforafrica.net/downloads/NRDS/Togo_fr.pdf), le Nigéria reste donc la seule destination possible de la ré-exportation du riz. D'ailleurs, presque l'entièreté des **"exportations officielles"** de riz au Bénin a pour unique destination le Nigéria. 

<div class="vizdiv" id="riceimport"> </div>

Les taxes douanières du Nigéria sur l'importation du riz blanchi [sont effectivement exhorbitantes](https://www.unifiedjournals.org/ujafs/pdf/2016/jane.pdf). Depuis quelques années, le pays maintient également une interdiction sur l'importation via les frontières terrestres et l'achat à partir de devises étrangères<sup> [3](https://www.export.gov/article?id=Nigeria-Prohibited-and-Restricted-Imports), [4](https://www.premiumtimesng.com/news/headlines/254114-nigeria-ban-rice-importation-year-buhari.html)</sup>. L'objectif visé ? atteindre une autosuffisance alimentaire en stimulant l'agriculture locale. Cette politique a permis une réduction phénoménale de l'importation du riz étuvé dès 2016<sup> [5](https://newtelegraphonline.com/2018/01/2017-rice-import-thailand-drops-23192mt/), [6](https://guardian.ng/news/nigerias-rice-importation-drops-to-20000mt/)</sup>. Bien que l'objectif initial, en ce qui concerne la production locale, ne soit pas encore atteint, une nette augmentation a quand même été observée<sup> [7](http://www.punchng.com/rice-production-in-nigeria-increases-to-5-8m-tonnes-in-2017-rifan/), [8](https://www.thisdaylive.com/index.php/2017/12/13/cbn-local-rice-production-records-70-increase/)</sup>. Il faut dire que le Nigeria a compris quelque chose en ce qui concerne l’autosuffisance alimentaire et pas nous...

<div class="midtitle">"un business très juteux"</div>

La politique nigériane sur l'importation du riz décortiqué fait [des heureux, au Cameroun et au Bénin](https://www.agenceecofin.com/riz/0304-18937-la-taxe-de-110-sur-les-importations-de-riz-fait-la-fortune-des-contrebandiers-beninois-et-camerounais). Un nouveau **secteur économique** très lucratif [a vu le jour](https://www.bloomberg.com/news/articles/2018-03-21/smugglers-run-riot-as-nigeria-tries-to-keep-foreign-rice-at-bay): `la contrebande de riz`. Le système est simple: tirer profit de la [Taxe Extérieur Commune (TEC), adoptée par le Bénin](http://douanes-benin.net/index.php/2017/09/26/tec-cedeao-sh-2017/) en tant que membre de la CEDEAO/UEMOA. Cette TEC fixe les droits de douane sur l'importation du riz à 10% (l'une des plus basse au monde), négligeable face au ~60% du Nigéria. Pour les contrebandiers, il suffit donc d'importer le riz via le port autonome de **Cotonou** puis faire passer clandestinement la cargaison à travers les frontières porreuses et corrompues qui nous séparent du Nigéria. Encore plus astucieux ? déclarer l'importation comme destinée au pays de l'hinterland pour être exempté de certaines taxes sur la consommation.

On pourrait tenter de dédramatiser en avançant que le trafic informel du riz occupe et nourrit la population désoeuvrée. Personnellement, je trouve cet argument fallacieux. Qu'on ne se voile pas la face ! La contrebande du riz est bien organisée. Elle ne profite par à la ménagère désespérée, mais plutôt à des businessmen véreux qui cherchent à duper le *système*.

<div class="midtitle">"tout pour entraver l'autosuffisance en riz"</div>

Même en ignorant *l'aspect contrebande*, l'importation massive de riz blanchi à d'autres répercussions toutes aussi alarmantes. En effet, l'importation destinée à la consommation, elle aussi en hausse depuis quelques années, [heurte notre secteur agricole déjà sous-developpée](https://www.agenceecofin.com/riz/1007-21466-cedeao-la-taxe-de-10-sur-le-riz-importe-ne-protegera-pas-la-production-locale) et interfère avec nos espoirs pour une autosuffisance alimentaire.

Étant donnée l'incapacité de la production locale à satisfaire, actuellement, la demande croissante de riz sur le territoire, il est vrai que l'importation est inévitable. Elle pourrait même rendre la riziculture locale plus compétitive. Cependant, maintenir des conditions qui ne contribuent qu'à inonder le marché de produits importés a sûrement l'effet contraire... La solution nigériane semble toutefois inadaptée et trop radicale, d'autant plus que nous avons d'autres problèmes à résoudre avant.

<div class="vizdiv" id="ricebenin"> </div>

Le Bénin dispose de conditions édaphiques et climatiques bien réparties sur l'ensemble du territoire et [très propices à la riziculture](http://www.inter-reseaux.org/IMG/pdf/Rapport_Final_Etat_des_Lieux_Riz_1_.pdf). D'après [les rapports du Conseil de Concertation des Riziculteurs du Benin (CCR-B)](http://www.inter-reseaux.org/IMG/pdf/Rapport_Final_Etat_des_Lieux_Riz_1_.pdf), la production locale aurait doublée entre 2007 et 2014, et les efforts menés pour améliorer le rendement et la qualité du produit final semblent également fructueux. Mais tout comme les autres secteurs de l'agriculture vivrière, la culture du riz est confrontée à un nouveau défit : [l'accaparement massif des terres agricoles exploitables par des sociétés étrangères, parfois en connivence avec certaines ONGs](https://www.pambazuka.org/fr/governance/des-paysans-béninois-disent-non-à-l’accaparement-des-terres). Dès acquisition, les cultures vivrières sur ces terres sont souvent remplacées par des projets de production de biocarburant<sup> [11](http://www.fondation-farm.org/zoe/doc/foncier_benin.pdf), [12](https://dumas.ccsd.cnrs.fr/dumas-00948184/document), [13](http://terres.redtac.org/IMG/pdf/revue_de_la_litterature_finale_-_version_i_-_pre-atelier-2.pdf)</sup> ([ex: culture de jatropha](https://www.wsj.com/articles/SB118788662080906716?mod=googlenews_wsj))


<img class="img-center" src="https://i.imgflip.com/28qlfp.jpg"/>

Un autre problème auquel doivent faire face les producteurs, est la réticence de la population à consommer le riz blanc local. Nos palais devenus bourgeois, préfèrent le *parfum exquis* du riz thaï/basmati asiatique à celui local, [qualifié à tort d'être de moins bonne qualité](https://youtu.be/OYKQMIrOf-Q?t=1m54s). 

Tant que des mesures appropriées ne seront pas prises, l'avenir du riz blanc béninois sera *tout le contraire de sa couleur*.


### Nécéssité d'un dialogue impliquant tous les acteurs. 

L'importance d'une autosuffisance en riz, l'État béninois en est bien conscient. Encore récemment, le secteur riz était l'une des priorités du gouvernement. Appuyé par le financement de l'Union Européenne à hauteur de 3.5 milliard de FCFA, entre 2009 et 2013<sup> [14](http://www.cantool.net/download/250/pafiriz.pdf), [15](http://www.rfi.fr/afrique/20170210-tres-chemin-croix-riz-beninois-local-importations-delice-nigeria)</sup>, le Bénin a fortement encouragé la filière. Comment comprendre alors pourquoi le gouvernement demeure sourd aux demandes parfaitement raisonnables des producteurs<sup> [16](http://www.rfi.fr/afrique/20170210-tres-chemin-croix-riz-beninois-local-importations-delice-nigeria), [17](https://westafrica.rikolto.org/en/node/1317)</sup>, en maintenant notamment une politique d'importation absurde ? 

<div class="midtitle">"nouvelle plateforme de promotion par la plateforme de promotion du riz ?"</div>

Le problème ne se limite pas seulement à l'absence de réformes au niveau de l'État, il faudrait aussi tenir responsable le conseil de concertation des riziculteurs qui a échoué à promouvoir efficacement ses produits. Ce n'est évidemment pas par manque d'avoir essayé, mais les démarches entreprises dans la direction ne peuvent être qualifiées que de `«circle jerking bureaucratique»`. 

En mai 2017, le CCR-B, en partenariat avec plusieurs professionnels des médias, avait constitué la [*"Plateforme d'informations et de Promotion du Riz du Benin (PiPr RB)"*](https://matinlibre.com/index.php/societe/item/11379-promotion-de-la-consommation-du-riz-made-in-benin-la-plateforme-des-professionnels-des-medias-constituee) ayant pour promesse: `«Œuvrer pour une meilleure adoption du riz béninois par les béninois, et leur Etat»`. Un an plus tard, on se demande bien où sont passés les résultats tant attendus de cette plateforme... Un véritable dialogue entre consommateurs, producteurs et autorités s'avère donc indispensable pour espérer un salut.

<div class="midtitle">"des solutions potentielles ignorées"</div>

C'est bien évident que nous avons besoin de nouvelles réformes facilitant l'accès à la production locale, et une règlementation plus réfléchie en ce qui concerne l'importation. Un simple relèvement de la taxe d'importation du riz à 20% pourrait en effet être [avantageux sur tous les plans](http://www.ecoasso.org/articles/Agossadou_et_al.pdf). 
Par ailleurs, puisque [notre économie dépends fortement du Nigeria](https://www.reuters.com/article/us-nigeria-benin-smuggling/nigeria-recession-deals-blow-to-smuggling-hub-benin-idUSKBN17125X) et que le déficit en production de riz de ce dernier demeure énorme, l'adoption de mesures conjointes sur la gestion du riz pourrait être bénéfique aux deux pays. Le Nigeria y gagnera une réduction de la contrebande de riz, et le Bénin, un partenaire commercial pour `l'exportation formelle` de la production locale. 

En attendant les réformes de l'État, nous, en tant que citoyens, pouvons déjà mener certaines actions simples mais significatives comme *consommer local*. Après tout, le développement de notre pays est une responsabilité partagée. 


<br>

## Sources

### Données et viz

[https://github.com/maclandrol/rice_benin](https://github.com/maclandrol/rice_benin)

### Selection de documentaires et lectures que je recommande.

- [Mini-documentaire: Riz du Bénin, riz de demain](https://www.youtube.com/watch?v=zCNCoCh8FjQ)
- [Reportage France 24 sur l'accaparement de terres arables](https://www.youtube.com/watch?v=GxFTGq94dXs)
- [Le secteur informel en Afrique de l'Ouest francophone](https://openknowledge.worldbank.org/bitstream/handle/10986/9364/9782744076602.pdf)
- [Discussion sur la ré-exportation vers le Nigéria, pilier de l'économie béninoise](http://horizon.documentation.ird.fr/exl-doc/pleins_textes/pleins_textes_7/autrepart/010014755.pdf)
- [Un article qui résume la situation du riz au Nigéria (en)](https://www.pmnewsnigeria.com/2018/03/13/the-rice-war-how-asian-rice-importers-sabotage-nigerias-rice-policy/)
- [Un autre qui parle de la contrebande (en)](https://www.bloomberg.com/news/articles/2018-03-21/smugglers-run-riot-as-nigeria-tries-to-keep-foreign-rice-at-bay)
- [Enquête complète publiée sur RFI: Le chemin de croix du riz béninois](http://www.rfi.fr/afrique/20170210-tres-chemin-croix-riz-beninois-local-importations-delice-nigeria)
- [Étude sur la compétitivité de la riziculture béninoise](http://www.hubrural.org/IMG/pdf/etude_competitivite_riziculture_beninoise.pdf)
- [CCR-B: État des lieux de la filière riz au Bénin en 2014](http://www.inter-reseaux.org/IMG/pdf/Rapport_Final_Etat_des_Lieux_Riz_1_.pdf)



<script src="/assets/data/rice/trade.js" type="text/javascript" charset="utf-8" ></script>
<script src='//d3plus.org/js/d3.js' type='text/javascript'></script>
<script src='//d3plus.org/js/topojson.js' type='text/javascript'></script>
<script src='//d3plus.org/js/d3plus.js' type='text/javascript'></script>
        
<script>
(function(){

	var mapcolor = d3.scale.ordinal()
    .domain([1, 95])
    .range(["#f0a40f", "#d20aee", "#5bbfe9", "#3bcd4d", "#ef1b64", "#a37060", "#5e70ed", "#e93903", "#86c194", "#ee97d6", "#798217", "#dd2ab1", "#777a9b", "#bc651f", "#d6ad6a", "#aaacff", "#258987", "#a2bf32", "#129031", "#c95189", "#c3afaa", "#b256c4", "#d54e4f", "#f59d81", "#1582c3", "#966bb0", "#847d4b", "#2ec7be", "#288c5e", "#7dc567", "#965bed", "#aeba68", "#fd95ac", "#88bdbf", "#9d761b", "#cdb301", "#33ca93", "#e8228a", "#777e73", "#a86988", "#c6a9d5", "#f4163d", "#fc9c54", "#6ec92f", "#247ed8", "#f88deb", "#d94c23", "#c35c62", "#bdb47f", "#68baff", "#dea696", "#7075c4", "#c7b351", "#59884a", "#95b7d4", "#dd3e76", "#d0a2ea", "#b949d9", "#d3419d", "#ddac37", "#b7664d", "#b25f9c", "#aa6f37", "#69835f", "#3a859b", "#a0b2ea", "#49c3d4", "#4b8b13", "#9065d8", "#c95b39", "#e39fc0", "#e63a3c", "#9e63c4", "#48cd2e", "#ff9b3b", "#f61209", "#9dbf4f", "#cb33d9", "#67c94d", "#2b79ec", "#eba553", "#5f8473", "#95bcaa", "#ae6875", "#897d35", "#ec1e77", "#d1ad80", "#608733", "#cc44b1", "#8b7874", "#e13c63", "#73c57d", "#967188", "#aeb6aa", "#96774c", "#52c6a9", "#bab0bf", "#5d809b", "#a8bb7f", "#687baf", "#8174af", "#b36761", "#7f7e5f", "#e8a56a", "#70834a", "#8db3ff", "#cb5a21", "#a3b6bf", "#ce5076", "#c3b469", "#e598eb", "#eb9eab", "#5877d8", "#5cc966", "#c3539d", "#daac52", "#b7b934", "#c65b4e", "#f19e96", "#517cc3", "#f99c6b", "#896dc4", "#a961b0", "#e43b50", "#cab335", "#d8a7ab", "#4781af", "#f696c1", "#afb1d5", "#ac4ded", "#d93f89", "#8d729c", "#d0a8c0"]);

	var viz_d3viz_87 = d3plus.viz()
	    .container('#tradeplot')
	    .type('line')
	    .color('Trade Flow')
		.text('Trade Flow')
		.legend({'font': {'size': 12},'data': false,'size': [20,50]})
		.y({'grid': false,'value': 'Trade Value (US$)'})
		.x({'ticks': {'labels': [1998,2000,2002,2004,2006,2008,2010,2012,2014,2016]},'grid': false,'value': 'Year'})
		.id('Trade Flow')
	    .data(trade_data)
	    .title({'value':"Importations et Exportations déclarée par le Bénin entre 1998 et 2016", 'sub': "Source: ONU/Comtrade"})
	    .draw();

    var viz_d3viz_94 = d3plus.viz()
        .container('#trademap')
        .data(import_dt_2016)
        .type("tree_map")
        .id(['L1','L2'])
        .size('Trade Value (US$)')
		.title({"value":"Biens et services importés en 2016 par le Bénin", 'sub': "Source: ONU/Comtrade"})
        .legend(false)
        .color({
	      "scale": mapcolor,
	      "value": "L1"
	    })
        .tooltip({'children': false})
        .depth(0)
		.footer({'value': '35.7 % des importations proviennent de la Thaïlande, l\'Inde et la Chine, essentiellement en raison de l\'importation massive de riz.', 'padding': 20, 'font': {'familly': '"Helvetica Neue", Lato, Arial, sans-serif;', 'size':12}})
		.draw();


	viz_d3viz_94.format({
	      "text": function(text, params) {
	        
	        if (text === "Trade Value (US$)") {
	          return "Trade Value";
	        }
	        else {
	          return d3plus.string.title(text, params);
	        }
	        
	      },
	      "number": function(number, params) {
	        
	        var formatted = d3plus.number.format(number, params);
	        
	        if (params.key === "Trade Value (US$)") {
	          return "$" + formatted + " USD";
	        }
	        else {
	          return formatted;
	        }
	     }})
        .ui([
        	{
        	"method" : function(value){
        		if (value == 'Import'){
					viz_d3viz_94.data(import_dt_2016)
					.title({"value":"Biens et services importés en 2016 par le Bénin", 'sub': "Source: ONU/Comtrade"})
					.footer({'value': '35.7 % des importations proviennent de la Thaïlande, l\'Inde et la Chine, essentiellement en raison de l\'importation massive de riz.', 'padding': 20, 'font': {'familly': '"Helvetica Neue", Lato, Arial, sans-serif;', 'size':12}});
        		}
        		else{
					viz_d3viz_94.data(export_dt_2016)
					.title({"value":"Biens et services exportés en 2016 par le Bénin",  'sub': "Source: ONU/Comtrade"})
					.footer({'value': '59,2 % des recettes d’exportations proviennent du coton brut et de l\'anacarde, principalement destinés à des pays de l\'Asie (Inde, Malaisie, Bangladesh, Chine)', 'padding': 20, 'font': {'familly': '"Helvetica Neue", Lato, Arial, sans-serif;', 'size':12}});

        		}
        		viz_d3viz_94.draw();
			},
        	"label": "Échanges",
        	"type": "toggle",
        	"font":{"size":15, "weight": 400},
        	"value"  : ["Import", "Export"]
      	}
	    ]);



    var viz_d3viz_2 = d3plus.viz()
        .container('#geodata')
        .data(geo_data)
        .type("geo_map")
        .coords({'value': '/assets/data/rice/countries.json', 'threshold': 0.2})
        .id('Country')
        .legend({'gradient':{"height": 8}, 'value':true, 'opacity': 1.0})
        .color({'heatmap':['#f0f9e8','#ccebc5','#a8ddb5','#7bccc4','#43a2ca','#0868ac'], 'value':'Riz importé ($US)'})
        .text('Pays')
        .title('Comparaison de la quantité de riz importée en 2016 dans une sélection de pays.')
        .tooltip(['Pays','Population','Riz importé ($US)','Quantité (kg)'])
        .footer({'value': 'Le Bénin importe bien plus de riz "per capita" que la plupart des autres pays. La carte montre, à titre comparatif le montant total des importations d\'une sélection de pays. Les pays de l\'Asie (souvent producteurs de riz) sont exclus', 'padding': 20, 'font': {'familly': '"Helvetica Neue", Lato, Arial, sans-serif;', 'size':12}})        
        .draw();


    var viz_rice_import = d3plus.viz()
        .container('#riceimport')
        .type('bar')
        .color('Reporter')
		.text('Reporter')
		.tooltip(['Trade Value (US$)','Netweight (kg)'])
		.legend({'font': {'size': 12},'data': false,'size': [20,40]})
		.y({'grid': false, 'value':'Trade Value (US$)'})
		.x({'grid': false, 'value':'Year'})
		.id('Reporter')
		.title({"value":'Importation de riz au Bénin, Togo et Nigeria', 'sub': "Source: ONU/Comtrade"})
        .data(rice_import)
        .footer({'value': 'L\'importation du riz au Togo reste constante tandis qu\'au Bénin, elle est en pleine croissance et connaît des pics à chaque fois qu\'un déclin est observé au Nigéria. Les importations du Nigéria en 2015 sont manquantes', 'padding': 20, 'font': {'familly': '"Helvetica Neue", Lato, Arial, sans-serif;', 'size':12}})
        .draw();

     var viz_ricebenin = d3plus.viz()
        .container('#ricebenin')
        .type('line')
        .color({'value':'Attribute_Description', 'scale':['#f77f6c', '#736ad8','#2DC48A']})
		.text('Attribute_Description')
		.legend({'font': {'size': 12},'data': false})		
		.y({'value': 'Total en tonne', 'grid':false})
		.x({'grid': false,'value': 'Année','ticks': {'labels': [1990,1995,2000,2005,2010,2015]}})
		.id('Attribute_Description')
		.title({"value":"Quantité de riz produite, importée et nécessaire au Bénin","sub":"Source: USDA/ CCR-B/ INSAE"})
        .footer({'value': 'Les importations de riz comblent largement le déficit engendré par la production insuffisante. Notez qu\'il s\'agit d\'une estimation très conservative des besoins, basée sur une consommation annuelle de 35 kg/an/habitant et sur des données historiques du CCR-B (voir code). Les importations ne concernent que celles destinées à la consommation locale (corrigées pour tenir compte de la ré-exportation informelle)', 'padding': 20, 'font': {'familly': '"Helvetica Neue", Lato, Arial, sans-serif;', 'size':12}})
		.data(benin_rice);

	viz_ricebenin.format({
	      "text": function(text, params) {
	        
	        if (text === "Total en tonne") {
	          return "Quantité (en tonne)";
	        }
	        else {
	          return text;
	        }
	        
	      },
	      "number": function(number, params) {
	        
	        var formatted = d3plus.number.format(number, params);
	        
	        if (params.key === "Total en tonne") {
	        	return formatted;
	        }
	        else {
	          return number;
	        }
	     }})
        .draw();

})();


</script>

