# cantilever-new

`modal3`와 `uniform3` 형상을 pygalmesh volume mesh → ANSYS SOLID187 TET10으로 변환하고, 세 재료의 무구속 free-free modal 해석을 실행합니다.

## 입력 형상

GitHub Release의 다음 파일을 내려받아 프로젝트 루트에 압축 해제합니다.

- `modal_cantilever3.obj.zip`
- `uniform_cantilever3.obj.zip`

보정된 최신 원본 3MF는 `geometry/`에 있습니다. Release의 OBJ는 현재 폴더에 있던 보정 전 변환본이므로, 최신 반지름을 해석하려면 `geometry/`의 3MF를 OBJ로 다시 변환해 교체해야 합니다. OBJ는 GitHub의 파일 크기 제한 때문에 저장소가 아닌 Release asset으로 배포합니다.

## 실행

Python 환경에 `numpy`, `meshio`, `pygalmesh`가 필요합니다.

```bash
python run.py
```

ANSYS를 자동으로 찾지 못하면:

```bash
python run.py --exe "C:\\Program Files\\ANSYS Inc\\v252\\ansys\\bin\\winx64\\ANSYS252.exe"
```

`run.py`가 다음 작업을 순서대로 수행합니다.

1. OBJ를 pygalmesh TET4 volume mesh로 변환
2. edge 중간 절점을 추가해 SOLID187 TET10 CDB 생성
3. `runs/`에서 Inconel 718, Structural Steel, Formlabs Nylon 12 해석
4. `runs/free_free_modal_frequencies_all.csv` 저장

경계조건은 없는 true free-free입니다. 추출된 12개 모드 중 1–6차는 강체 모드, 7–12차는 탄성 고유진동 `f1–f6`입니다.
