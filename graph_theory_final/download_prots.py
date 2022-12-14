import wget
cat_name = "NEUROTRANSMITTER"

# ids = "P05177,P10635,P10632,Q96C23,Q9HCG7,P04062,P47712,P08684,P05108,O15269,O95237,O95470,P07902,P22570,P49619,Q5VZY2,Q6ZNC8,Q9HAY6,Q9HBH5,Q9NUN7,Q02928,P11712,Q6GTS8,O14495,P33260,P51589,P52429,Q16760,Q53H12,Q5KSL6,Q8NEB5,Q8NFR3,Q8NFU5,Q8TDN7,Q92781,Q969W0,Q9HA82,Q9NUV7,Q9Y6T7,Q9NZ01,Q16678,A5PLL7,Q02083,Q8N9I5,Q96SQ9,Q9UJ83,O60218,O75912,P52824,Q13574,Q14376,Q16739,Q643R3,Q6P531,Q7Z449,Q8N5B7,Q8NBN7,Q9HCS2,Q9UHH9,Q9UJ14,O14494,O75907,P00352,P47895,P50053,P51687,Q06136,Q6P1A2,Q6ZMG9,Q6ZWT7,Q8NF37,Q96G23,Q96N66,Q96NR8,Q9NY59,Q9UHE5,Q9UHK6,P04798,A6NGU5,O15270,O43688,O60906,O75452,O94788,P23743,P27544,P49675,P51570,Q13510,Q5QJU3,Q7L5N7,Q86XP1,Q8IU89,Q8IZV5,Q8N3Y7,Q8TC12,Q99999,Q9NR71,O00154,P98187".split(",")
# ids = "Q9NY46,Q9Y5Y9,Q9UI33,Q15858,O00180,Q7Z418,Q9Z2T2,Q91WD2,Q4KMQ2,Q7Z3S7,P61981,Q6UXB4,P50148,Q2PKF4,Q12965,O00160,Q8WXR4,Q9Y2K3,Q9Y623,Q9Y4I1,Q13459,Q9HD67".split(",")

# ION_CHANNELS
# https://www.uniprot.org/uniprotkb?query="ion%20channel"%20AND%20%28length%3A%5B200%20TO%201000%5D%29
# ids = "Q7NDN8,Q401N2,Q5H8A6,Q9LNJ0,Q9SKD6,Q9LEQ3,Q8GWD2,O65717,O82226,Q8RWS9,Q94AS9,O00299,Q09917,Q95Y52,Q9U358,Q9Y696,Q8C1E7,O65718,Q6NXK8,Q708S7,Q708S6,Q708S8,Q9Z0W7,Q9Z1Q5,Q9XSA7,Q9BXJ8,Q84W41,Q8LGN1,O14050,Q96FT7,Q9M8W7,P78348,P06971,Q24278,P46098,Q8WWG9,Q6X1Y6,G5ECT0,O35240,Q19351,Q1XA76,Q5H8A5,Q9UHC3,Q9C8E7,P55926,Q93YT1,Q8GXJ4,Q56X46,Q8L7Z0,P39719,Q5N941,Q9SJA4,Q9SL29,Q9SU64,Q9LD40,Q9M0A4,Q9S9N5,P35563,P23979,Q9SKD7,Q866Y9,P54245,Q9FXH6,Q16515,Q925H0,A8XNX8,O70212,Q708S3,Q708S4,Q4VY51,Q62962,O81078,Q9FH75,Q9LFN5,Q9SHV1,Q9SHV2,Q8LGN0,Q9LFN8,Q9LV72,O81776,Q6RHR6,Q9C5V5,Q08967,Q9NY37,A0A072VMJ3,G7IBJ4,G7JND3,Q9LDR2,Q7XP59,Q708S5,Q9LD37,Q7XJL2,Q7T1N4,Q9SW97,O04660,Q9SDQ4,Q9JHS6,F4IME2,Q84M97,O60741,P07510,P14867,P17787,P23416,P43681,P47870,Q05586,Q16281,Q8WXA8,Q9NQW8,O95264,P30926,P36544,P39086,P42263,P48058,Q04844,Q13002,Q15822,Q16280,Q70Z44,Q99572,Q9UGM1,P02708,P23415".split(",")

# ENZYMES
# https://www.uniprot.org/uniprotkb?query=enzyme%20AND%20(organism_id:9606)%20AND%20(length:[200%20TO%201000])
# ids = "P23368,P49427,Q712K3,Q9H0E7,P62068,Q70CQ1,Q70EK9,Q86T82,Q04446,Q8IU60,Q8TBC4,Q9P109,O00303,P48163,Q9GZZ9,Q9UK59,Q16798,Q9BYF1,Q8N2K1,Q8IZD4,P61086,Q9UBE0,Q9Y385,Q8IV48,Q9UBT2,P0DPD6,P09958,P42892,O95352,Q16763,Q9H832,Q5VVX9,Q969T4,Q8NBK3,Q8IUX4,Q8WVN8,Q08426,P41238,Q96LR5,P0DPD8,P48147,Q8NBJ7,Q9BWT3,P16083,P49768,Q86VQ3,Q8N608,Q8NET6,Q8WVM0,Q8WWY8,Q9BQ52,Q9H5Q4,Q9NUW8,Q9Y6K0,Q9Y6Y8,O43529,P15374,P27695,Q70EL3,Q86SR1,Q8IXK2,Q8NCW6,Q96EB6,Q96JF0,Q96JJ7,Q9Y646,O95395,O95396,O95453,Q5K4E3,Q5MY95,Q5VVQ6,Q7L1S5,Q8IVS2,Q8NFW8,Q9BY49,Q9H2A9,Q9H777,Q9H8X2,Q9NPH2,Q9Y251,O43272,P30793,Q66K79,Q6P179,Q6PIY7,Q6Y288,Q8NBQ5,Q8NCH0,Q8NCL4,Q8TCT0,Q8WUD6,Q92542,Q92560,Q96BI3,Q96D53,Q96GR2,Q9BRJ7,Q9NQZ7,Q9NRB3".split(",")

# PROTEIN KINASES, apparnently...
# https://www.uniprot.org/uniprotkb?query=(cc_bpcp_temp_dependence:5)
# ids = "Q15067,P00505,P45563,Q8T1G4,Q84W55,O32224,P0ADG4,P26364,P0A9B6,P39065,Q12068,Q8VIJ5,Q9V4C8,P11759,P9WFZ5,Q41131,P80561,D4GYV5,F1C7I4,Q17CS8,C0HJB3,Q9HVW0,Q75UV1,Q6TGQ8,Q1GNW5,B5BLW5,P17556,Q5L2C2,Q2TDY4,L8B068,Q59196,G4FEF4,A4ECA9,Q3L8N0,Q9AQS0,O28007,R4ZGN4,F8G0M4,Q9WYH1,K9NBS6,D4GUJ7,J7GQ11,I2DBY1,Q4PNI0,Q45515,Q9WYR4,B4XC07,Q57679,Q58487,Q65CX4,Q9RFR0,Q8NK92,O28608,Q9KI47,Q868M7,G9BY57,Q45070,T2KNC2,O75003,Q8RJB2,A0A1S6EK91,O57936,P27142,E1VA04,Q9K4Y9,Q76KX0,B3EWG2,P9WIT3,Q6KZ25,B7GPC7,O85078,O28523,Q9LCU3,P80305,Q8DK23,O86959,K7N5L0,P85513,C5B2R8,F5L9G2,P85147,Q9Z8L0,P49421,P17557,Q4WZ11,B9LW38,P86833,Q5FJ41,P21800,P82288"

# NEUROTRANSMITTER
# https://www.uniprot.org/uniprotkb?query=neurotransmitter%20AND%20(length:[0%20TO%201000])
ids = "O76689,Q9H1V8,Q8BJI1,P31662,Q9H2J7,G5EBN9,Q08469,Q8BG16,Q9XS59,Q5R9C2,Q9GZN6,Q28001,O14804,P22303,P31645,P48066,P61266,P80404,Q05329,Q13277,Q16568,Q16572,Q496J9,Q9UJD0,Q9Y215,P48065,Q12846,Q16623,Q7Z7G2,Q8WVH0,Q8WZ04,Q99259,Q9P2U7,Q9P2U8,Q9Y345,P21397,Q05084,Q7L0J3,Q8NDX2,Q96N87,Q99884,Q9H598,O14810,P07101,P20366,P23975,P28329,P31641,P54219,Q05940,Q7L1I2,Q9NSD5,P01213,P21964,P30531,P48067,Q01959,Q13519,Q6PUV4,Q9GZV3,G5EBM5,A0A193KX02,O16000,O35167,O35304,O35526,O55192,P21398,P23978,P24529,P31646,P38433,P50554,P63040,P90986,P97411,Q02563,Q62923,Q63054,Q761V0,Q810F0,Q86B61,Q9JI12,Q9JIS5,O35417,O97148,P18088,P21396,P21836,P28571,P28573,P31650,P58295,P61264,Q3TXX4,Q62634,Q7KVY7,Q8BG39,Q9JIR3,Q9MZ34".split(",")


for ID in ids:
    try:
        wget.download(
            f"https://alphafold.ebi.ac.uk/files/AF-{ID}-F1-model_v2.cif",
            f"./inp/{cat_name}/AF-{ID}.cif")
    except:
        print("no good!", ID)

