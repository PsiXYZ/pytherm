import pytherm.activity.pitzer as pz


def test_get_y_harvie84():
    t_rez = {
        'Na': 0.8416, 'Cl': 0.7930, 'OH': 0.8932
    }
    m_nacl = 2
    m_naoh = 2
    ph = {
        "Na": m_nacl + m_naoh,
        "Cl": m_nacl,
        "OH": m_naoh,
    }
    am = pz.Pitzer(ph, pz.datasets.Harvie84)

    y = am.get_y(ph)
    for i in y:
        if round(y[i], 4) == t_rez[i]:
            continue
        else:
            assert False
    assert True
