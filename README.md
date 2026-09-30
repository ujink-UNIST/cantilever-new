# cantilever-new

`modal3`와 `uniform3` 형상을 pygalmesh volume mesh → ANSYS SOLID187 TET10으로 변환하고, 세 재료의 무구속 free-free modal 해석을 실행합니다.

## 입력 형상

OBJ는 용량 때문에 Git 저장소가 아닌 `v1.0` Release에 압축해 두었습니다. `run.py`는 OBJ가 없으면 GitHub CLI로 다음 파일을 자동 다운로드·압축 해제합니다.

- `modal_cantilever3.obj.zip`
- `uniform_cantilever3.obj.zip`

private 저장소이므로 워크스테이션에서 한 번 `gh auth login`이 필요합니다. 최종 OBJ 실체 부피는 modal `83,319.77 mm³`, uniform `83,341.32 mm³`로 차이는 `0.026%`입니다. uniform 공통 반지름은 `0.806104 mm`입니다.

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
